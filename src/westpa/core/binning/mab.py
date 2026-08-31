import logging
from os.path import expandvars
from typing import List, Optional

import numpy as np

import westpa
from westpa.core.binning import FuncBinMapper
from westpa.core.binning.assign import index_dtype, rectilinear_assign_python

log = logging.getLogger(__name__)


class MABBinMapper(FuncBinMapper):
    """
    Adaptively place bins between minimum and maximum segments along
    the progress coordinate. Extrema and bottleneck segments are assigned
    to their own bins.
    """

    def __init__(
        self,
        nbins: List[int],
        direction: Optional[List[int]] = None,
        skip: Optional[List[int]] = None,
        bottleneck: int = 1,
        pca: bool = False,
        mab_log: bool = False,
        bin_log: bool = False,
        bin_log_path: str = "$WEST_SIM_ROOT/binbounds.log",
        strict_Z: bool = True,
    ):
        """
        Parameters
        ----------
        nbins : list of int
            List of number of bins in each dimension.
        direction : Optional[list of int], default: None
            List of directions in each dimension. Direction options:
                0   : default split at leading and lagging boundaries
                1   : split at leading boundary only
                -1  : split at lagging boundary only
                86  : no splitting at either leading or lagging boundary (both bottlenecks included)
        skip : Optional[list of int], default: None
            List of skip flags for each dimension. Default None (no skipping).
        bottleneck : int, default: 1
            Whether to enable bottleneck walker splitting. By default, one bottleneck segment will be chosen.
        pca : bool, default: False
            Whether to perform PCA on progress coordinates before bin assignment.
        mab_log : bool, default: False
            Whether to output MAB info to west.log.
        bin_log : bool, default: False
            Whether to output MAB bin boundaries to a log file.
        bin_log_path : str, default: "$WEST_SIM_ROOT/binbounds.log"
            Path to output bin boundaries.
        strict_Z : bool, default: True
            Whether to put bottleneck-like segments (highest Z value but not technically a bottleneck, i.e, Z < 0)
            into bottleneck bins or not.
        """
        # Verifying parameters
        if nbins is None:
            raise ValueError("nbins is missing")
        ndim = len(nbins)

        direction = direction or [0] * ndim
        if len(direction) != ndim:
            direction = [0] * ndim
            log.warning("Direction list is not the correct dimensions, setting to defaults.")

        skip = skip or [0] * ndim
        if len(skip) != ndim:
            skip = [0] * ndim
            log.warning("Skip list is not the correct dimensions, setting to defaults.")

        kwargs = dict(
            nbins_per_dim=nbins,
            direction=direction,
            skip=skip,
            bottleneck=bottleneck,
            pca=pca,
            mab_log=mab_log,
            bin_log=bin_log,
            bin_log_path=bin_log_path,
            strict_Z=strict_Z,
        )

        n_total_bins = self.determine_total_bins(**kwargs)

        super().__init__(map_mab, n_total_bins, kwargs=kwargs)

    def determine_total_bins(
        self, nbins_per_dim: List[int], direction: List[int], skip: List[int], bottleneck: int, **kwargs
    ) -> int:
        """
        Calculate the total number of bins needed, taking direction and skipping into account.
        This function is necessary because functional bin mappers need to "reserve"
        bins and tell the sim manager how many bins they will need to use, this is
        determined by taking all direction/skipping info into account.

        Parameters
        ----------
        nbins_per_dim : list of int
            Number of total bins in each dimension within the linear portion.
        direction : list of int
            Direction in each dimension.
        skip : list of int
            List indicating whether to skip each dimension.
        bottleneck : int
            Whether to include a separate bin for bottleneck walker(s).
        **kwargs : dict
            Additional MAB parameters (unused).

        Returns
        -------
        n_total_bins : int
            Number of total bins.

        """
        # Update nbins_per_dim with any skipped dimensions, setting number of bins along skipped dimensions to 1
        skip = np.array([bool(s) for s in skip])
        nbins_per_dim = np.array(nbins_per_dim)
        nbins_per_dim[skip] = 1

        # Total bins is product of all linear bins plus and special bins
        n_total_bins = nbins_per_dim.prod()
        for direct, skip_dim in zip(direction, skip):
            if not skip_dim:
                if direct in [-1, 1]:
                    # 1 lead or lag bin + 1 bottleneck bin
                    n_total_bins += 1 + 1 * bottleneck
                elif direct == 0:
                    # 2 lead/lag bins + 2 bottleneck bins
                    n_total_bins += 2 + 2 * bottleneck
                elif direct == 86:
                    # 0 lead/lag + 2 bottleneck bins
                    n_total_bins += 2 * bottleneck
        return n_total_bins


def map_mab(coords: np.ndarray, mask: np.ndarray, output: np.ndarray[index_dtype], *args, **kwargs) -> np.ndarray[index_dtype]:
    """
    Adaptively place bins based on extrema and bottleneck segments along the progress coordinate.

    Bottleneck segments are where the difference in probability is the greatest
    along the progress coordinate. Operates per dimension (unless skipped) and places a fixed number of
    evenly spaced bins between the segments with the min and max pcoord values. Extrema and
    bottleneck segments are assigned their own bins.

    Parameters
    ----------
    coords : np.ndarray
        An array with pcoord and weight info.
    mask : np.ndarray
        Boolean array to filter out unwanted segments.
    output : np.ndarray[index_dtype]
        The main array that, for each segment, holds the bin assignment.
    *args : list
        Additional arguments.
    **kwargs : dict
        Additional keyword arguments. Contains most of the MAB-needed parameters.

    Returns
    ------
    output : np.ndarray[index_dtype]
        Array with bin assignments for each segment.

    """
    # Argument Processing
    nbins_per_dim = kwargs.get("nbins_per_dim")
    ndim = len(nbins_per_dim)
    pca = kwargs.get("pca", False)
    bottleneck = kwargs.get("bottleneck", 1)
    direction = kwargs.get("direction", [0] * ndim)
    skip = kwargs.get("skip", [0] * ndim)
    mab_log = kwargs.get("mab_log", False)
    bin_log = kwargs.get("bin_log", False)
    bin_log_path = kwargs.get("bin_log_path", "$WEST_SIM_ROOT/binbounds.log")
    strict_Z = kwargs.get('strict_Z', True)
    binbounds_determination_mask = kwargs.get('binbounds_determination_mask', None)

    if not np.any(mask):
        return output

    if skip is None:
        skip = [0] * ndim

    allcoords = coords.copy()
    allmask = mask.copy()

    report = True if coords[-1, -1] == 1 else False  # Report only when binning final
    splitting = True if coords[-1, -1] == 1 else False  # Only split when binning final
    strict = True if binbounds_determination_mask is not None else False

    # Mask out everything not needed for bin boundary determination (min/max/bottleneck).
    # The "if" condition is a way to bypass the automatic behavior of using the same mask to both determine bin bounds + assign.
    # The `mask` is recommended to be a subset of `binbounds_determination_mask` but guard rails are removed (when latter is provided)
    # so points located outside of minlist/maxlist as determined by `binbounds_determination_mask` are clipped to the nearest bin.
    if strict:
        coords = allcoords[binbounds_determination_mask, :ndim]
        weights = allcoords[binbounds_determination_mask, ndim] if allcoords.shape[1] >= ndim else None
        mask = allmask[binbounds_determination_mask]
    else:
        coords = allcoords[allmask, :ndim]
        weights = allcoords[allmask, ndim] if allcoords.shape[1] >= ndim else None
        mask = allmask[mask]

    originalcoords = np.copy(coords)
    if pca and len(output) > 1:
        coords = apply_pca(coords, weights)

    # Computing special bins (bottleneck and boundary bins)
    minlist, maxlist, bottlenecks_forward, bottlenecks_reverse = calculate_bin_boundaries(
        originalcoords, weights, mask, skip, splitting, bottleneck, strict_Z
    )

    if mab_log and report:
        log_mab_stats(minlist, maxlist, direction, skip)

    # Assign segments to bins
    output, n_bottleneck_filled = bin_assignment(
        allcoords,
        allmask,
        minlist,
        maxlist,
        bottlenecks_forward,
        bottlenecks_reverse,
        nbins_per_dim,
        direction,
        skip,
        splitting,
        bottleneck,
        output,
        strict,
    )

    # Report MAB bin statistics
    if bin_log and report and westpa.rc.sim_manager.n_iter:
        log_bin_boundaries(
            skip,
            bottleneck,
            direction,
            bin_log_path,
            minlist,
            maxlist,
            nbins_per_dim,
            n_bottleneck_filled,
            bottlenecks_forward,
            bottlenecks_reverse,
            strict_Z,
        )

    return output


def apply_pca(coords, weights):
    colavg = np.mean(coords, axis=0)
    varcoords = coords - colavg
    covcoords = np.cov(varcoords.T, aweights=weights)
    eigval, eigvec = np.linalg.eigh(covcoords)
    eigvec = eigvec[:, np.argmax(np.abs(eigvec), axis=1)]
    eigvec[:, np.diag(eigvec) < 0] *= -1
    return np.dot(varcoords, eigvec)


def calculate_bin_boundaries(coords, weights, mask, skip, splitting, bottleneck, strict_Z):
    """
    This function calculates minima, maxima, and bottleneck segments.
    """
    skip = np.array([bool(s) for s in skip])

    # Initialize lists to hold bottleneck segments along each dimension
    bottlenecks_forward, bottlenecks_reverse = [None] * len(coords[0]), [None] * len(coords[0])
    # number of unmasked coords
    n_coords = mask.sum()
    # Grabbing all unmasked coords and weights
    unmasked_coords = coords[mask, :]
    unmasked_weights = weights[mask] if weights is not None else None
    # Replace any zero weights with non-zero values so that log(weight) is well-defined
    if unmasked_weights is not None:
        unmasked_weights[unmasked_weights == 0] = 10**-323

    # We calculate the min and max pcoord along each dimension (boundary segments) even if skipping
    maxlist = np.max(unmasked_coords, axis=0)
    minlist = np.min(unmasked_coords, axis=0)

    # Looping over each dimension of progress coordinate to calculate bottleneck
    for n in range(len(coords[0])):
        if splitting and bottleneck and not skip[n]:
            bottlenecks_forward[n], bottlenecks_reverse[n] = detect_bottlenecks(
                unmasked_coords, unmasked_weights, n_coords, n, bottleneck, strict_Z
            )

    # raise ValueError(f'{bottlenecks_forward=}, {bottlenecks_reverse=}')

    return minlist, maxlist, bottlenecks_forward, bottlenecks_reverse


def detect_bottlenecks(unmasked_coords, unmasked_weights, n_coords, n, n_bottlenecks, strict_Z):
    """
    Detect the bottleneck segments along the given coordinate n, this uses the weights
    """
    # Grabbing all unmasked coords in current dimension, plus corresponding weights
    # Sort by current dimension in coord, smallest to largest, then by weights
    sorted_indices = np.lexsort(
        (
            unmasked_weights,
            unmasked_coords[:, n],
        )
    )
    # sorted_indices = np.lexsort((unmasked_coords[:, n],))  # same as argsort code before

    # Grab sorted coords and weights
    coords_srt = unmasked_coords[sorted_indices, :]
    weights_srt = unmasked_weights[sorted_indices]

    # Short circuit out and return empty lists if only 1 or less segments
    # Those will be considered by the leading/trailing walkers
    if len(weights_srt) < 2:
        return [], []

    # Also sort in reverse order for opposite direction bottlenecks (reverse-stable), where equivalents
    # in the first key (n-th dimension coord) are already sorted by secondary key (weights).
    # Verified solution from https://stackoverflow.com/a/64243103
    flip_indices = (len(coords_srt) - 1) - np.argsort(coords_srt[::-1, n], kind='stable')[::-1]
    coords_srt_flip = coords_srt[flip_indices]
    weights_srt_flip = weights_srt[flip_indices]

    # Initialize the max directional differences along current dimension as None (these may not be updated)
    bottleneck_coords, bottleneck_coords_flip = None, None

    # Summing up weights of all walkers ahead of current walker along current dim in both directions
    # Starting from 1 by default because we care about what is ahead of the boundary walker (not including it).
    # Cumsum of the opposite direction starting from the second point, then reversing it.
    cumulative_prob = np.cumsum(weights_srt[1:][::-1])[::-1]
    cumulative_prob_flip = np.cumsum(weights_srt_flip[1:][::-1])[::-1]

    # Calculating the bottlneck walker based on difference of log weight of current walker
    # and cumulative weight of everything ahead (Z in the MAB paper).
    # We use the log as weights vary over many orders of magnitude.
    # Note a negative Z indicates the cumulative weight ahead of the current walker is larger than the weight of the current walker,
    # while a positive Z indicates the cumulative weight ahead of the current walker is smaller, indicating a barrier.
    # Efficiency is better with np.argmax when looking for one bottleneck walker.
    # Skips leading/trailing walker because Z is literally undefined for those points.
    if n_bottlenecks == 1:
        # Forward direction (coord -> +inf)
        # Calculate Z
        Z_array = np.log(weights_srt[:-1]) - np.log(cumulative_prob)
        Zmax_idx = np.argmax(Z_array)
        Zmax_value = Z_array[Zmax_idx]
        # If strict_Z, only pick segment as bottleneck if Z > 0 or skip, otherwise pass most bottleneck-like
        bottleneck_coords = [tuple(coords_srt[Zmax_idx, :])] if not strict_Z or Zmax_value > 0 else []

        # Do same for reverse direction (coord -> -inf)
        Z_array = np.log(weights_srt_flip[:-1]) - np.log(cumulative_prob_flip)
        Zmax_idx = np.argmax(Z_array)
        Zmax_value = Z_array[Zmax_idx]
        bottleneck_coords_flip = [tuple(coords_srt_flip[Zmax_idx, :])] if not strict_Z or Zmax_value > 0 else []
    elif n_bottlenecks > 1:
        # Stable sort (secondary index by weight) to query the n-largest weight in the forward direction.
        # Tries to get as many unique bins as possible, up to requested (n_botlenecks).
        Z = np.log(weights_srt[:-1]) - np.log(cumulative_prob)
        sorted_Z_idx = np.argsort(Z, kind='stable')[::-1]
        sorted_Z = Z[sorted_Z_idx]

        bottleneck_coords = set(
            [tuple(coords_srt[sorted_Z_idx[idx]]) for idx, Z in enumerate(sorted_Z[:n_bottlenecks]) if not strict_Z or Z > 0]
        )

        for bn_idx in range(n_bottlenecks, len(cumulative_prob)):
            if len(bottleneck_coords) == n_bottlenecks or (strict_Z and sorted_Z[bn_idx + 1] <= 0):
                # We have our list or no more Z > 0 left
                break
            elif not strict_Z or sorted_Z[bn_idx + 1] > 0:
                # Will break out cleanly even if running through the for loop through completion because of the following correspondance
                # len(cumulative_prob) + 1 == len(sorted_Z) + 1 == len(coords_srt) == len(weights_srt)
                # cumulative_prob[i] <=> sorted_Z[i] <==> coords_srt[i-1] <==> weights_srt[i-1]
                bottleneck_coords.add(tuple(coords_srt[sorted_Z_idx[bn_idx + 1]]))

        # Doing the same for the opposite direction
        Z = np.log(weights_srt_flip[:-1]) - np.log(cumulative_prob_flip)
        sorted_Z_idx = np.argsort(Z, kind='stable')[::-1]
        sorted_Z = Z[sorted_Z_idx]

        bottleneck_coords_flip = set(
            [
                tuple(coords_srt_flip[sorted_Z_idx[idx] + 1])
                for idx, Z in enumerate(sorted_Z[:n_bottlenecks])
                if not strict_Z or Z > 0
            ]
        )

        for bn_idx in range(n_bottlenecks, len(cumulative_prob_flip)):
            if len(bottleneck_coords_flip) == n_bottlenecks or (strict_Z and sorted_Z[bn_idx + 1] <= 0):
                break
            elif not strict_Z or sorted_Z[sorted_Z_idx[bn_idx + 1]] > 0:
                bottleneck_coords_flip.add(tuple(coords_srt_flip[sorted_Z_idx[bn_idx + 1]]))

    # Return sorted version, small to large because sets were unordered
    return sorted(bottleneck_coords), sorted(bottleneck_coords_flip)


def log_mab_stats(minlist, maxlist, direction, skip):
    with np.printoptions(legacy='1.25'):
        westpa.rc.pstatus("################ MAB stats ################")
        westpa.rc.pstatus(f"minima in each dimension:      {minlist}")
        westpa.rc.pstatus(f"maxima in each dimension:      {maxlist}")
        westpa.rc.pstatus(f"direction in each dimension:   {direction}")
        westpa.rc.pstatus(f"skip in each dimension:        {skip}")
        westpa.rc.pstatus("###########################################")
        westpa.rc.pflush()


def bin_assignment(
    coords,
    mask,
    minlist,
    maxlist,
    bottlenecks_forward,
    bottlenecks_reverse,
    nbins_per_dim,
    direction,
    skip,
    splitting,
    bottleneck,
    output,
    strict,
):
    """
    Assign segments to bins based on the minima, maxima, and
    bottleneck segments along the progress coordinate.
    """
    # Update nbins_per_dim with any skipped dimensions, setting number of bins along skipped dimensions to 1
    skip = np.array([bool(s) for s in skip])
    nbins_per_dim = np.array(nbins_per_dim)
    nbins_per_dim[skip] = 1
    direction = np.array(direction)
    ndim = len(nbins_per_dim)

    # Some pre-calculated stats for tracking occupied bottleneck bins
    nbn_forward = [len(bn) if bn is not None else 0 for bn in bottlenecks_forward]
    nbn_reverse = [len(bn) if bn is not None else 0 for bn in bottlenecks_reverse]
    n_bottleneck_filled = np.zeros((sum(nbn_forward) + sum(nbn_reverse)), dtype=bool)

    # Boolean arrays that track use of special bins along each dimension
    skip_bneck_fwd = np.array([d == -1 if bottleneck else True for d in direction]) + skip
    skip_bneck_rev = np.array([d == 1 if bottleneck else True for d in direction]) + skip
    skip_lead = np.array([d in [86, -1] for d in direction]) + skip
    skip_lag = np.array([d in [86, 1] for d in direction]) + skip

    # List of dimensions that are not skipped
    active_dims = np.array([n for n in range(ndim) if not skip[n]])

    # Compute the boundary bin ID offsets
    # In forward direction, this is all the linear bins
    boundary_bin_id_offset_fwd = nbins_per_dim.prod()
    # In reverse, we add the number of forward boundary bins to the offset
    boundary_bin_id_offset_rev = boundary_bin_id_offset_fwd + (~skip_lead).sum()

    # Compute the bottleneck bin ID offsets
    # In forward direction, bin IDs are offset by all linear and boundary bins
    bneck_bin_id_offset_fwd = boundary_bin_id_offset_rev + (~skip_lag).sum()
    # In reverse, we add the number of forward bottleneck bins to the offset
    bneck_bin_id_offset_rev = bneck_bin_id_offset_fwd + (~skip_bneck_fwd).sum() * bottleneck

    # Calculate the rectilinear bin bounds ahead of time.
    # Create small-width bins if minlist[i] == maxlist[i]
    bin_bounds = [
        (
            np.linspace(minlist[i], maxlist[i], nbins_per_dim[i] + 1)
            if minlist[i] != maxlist[i]
            else np.linspace(minlist[i], maxlist[i] + 0.1, nbins_per_dim[i] + 1)
        )
        for i in range(ndim)
    ]

    # Assign everything in linear bins first, all at once.
    # If binning final coords, then we don't really care about what is being assigned to the initial coordinates.
    output = rectilinear_assign_python(coords[:, :ndim], mask=mask, output=output, boundaries=bin_bounds, strict=strict)
    # rectilinear_assign(np.asarray([coords[i, :ndim]], dtype=np.float32), mask=np.asarray([mask[i]], dtype=bool), output=output, boundaries=bin_bounds, boundlens=bound_lens)
    # [bin_id] = temp_output

    # Loop through all walkers and overwrite bin id for specials (bottleneck or leading walker)
    for i in range(len(output)):
        # Skip masked walkers, these walkers bin IDs are unchanged
        if not mask[i]:
            continue
        # Initialize bin ID and special tracker for current coord
        # The special variable indicates a boundary or bottleneck walker (not assigned to the linear space)
        bin_id, special = -1, False

        # Searching for bottleneck bins first
        if splitting and bottleneck:
            for i_acdim, n in enumerate(active_dims):
                # Grab coord(s) of current walker
                coord = coords[i, :ndim]
                # Assign bottlenecks, taking directionality into account
                # Check both directions when using 0 or 86
                # Note: 86 implies no leading or lagging bins, but does add bottlenecks for *both* directions when bottleneck is enabled
                # Note: When strict_Z = False, all bottleneck bins will typically be filled unless a walker is simultaneously in bottleneck bins along multiple dimensions
                # or there are too few walkers to compute free energy barriers
                for bfid, bforward in enumerate(bottlenecks_forward[n]):
                    if (coord == bforward).all() and not skip_bneck_fwd[n]:
                        bin_id = bneck_bin_id_offset_fwd + (i_acdim * bottleneck) + bfid
                        n_bottleneck_filled[sum(nbn_forward[:i_acdim]) + bfid] = 1
                        special = True
                        break
                for brid, breverse in enumerate(bottlenecks_reverse[n]):
                    if (coord == breverse).all() and not skip_bneck_rev[n]:
                        bin_id = bneck_bin_id_offset_rev + (i_acdim * bottleneck) + brid
                        n_bottleneck_filled[sum(nbn_forward) + sum(nbn_reverse[:i_acdim]) + brid] = 1
                        special = True
                        break

        # Now check for boundary walkers, taking directionality into account
        # This should only be done after fully checking for bottleneck walkers
        if splitting and not special:
            for n in active_dims:
                # Grab coord of current walker along current dimension
                coord = coords[i, n]
                if (coord == maxlist[n]) and not skip_lead[n]:
                    bin_id = boundary_bin_id_offset_fwd + n - skip_lead[:n].sum()
                    special = True
                    break
                elif (coord == minlist[n]) and not skip_lag[n]:
                    bin_id = boundary_bin_id_offset_rev + n - skip_lag[:n].sum()
                    special = True
                    break

        # output is the main array that, for each segment, holds the bin assignment
        # Only rewrite if the bin_id actually changed (bin_id >= 0) and in a special bin.
        if special and bin_id >= 0:
            output[i] = bin_id

    return output, int(sum(n_bottleneck_filled))


def log_bin_boundaries(
    skip,
    bottleneck,
    direction,
    bin_log_path,
    minlist,
    maxlist,
    nbins_per_dim,
    n_bottleneck_filled,
    bottlenecks_forward,
    bottlenecks_reverse,
    strict_Z,
):
    ndim = len(nbins_per_dim)
    skip = np.array([bool(s) for s in skip])
    active_dims = np.array([n for n in range(ndim) if not skip[n]])
    max_bottleneck = np.sum([bottleneck if direction[n] in [-1, 1] else 2 * bottleneck for n in active_dims]) if bottleneck else 0
    with open(expandvars(bin_log_path), 'a') as bb_file, np.printoptions(legacy='1.25'):
        # Iteration Number
        bb_file.write(f'Iteration: {westpa.rc.sim_manager.n_iter}\n')
        bb_file.write('MAB linear bin boundaries: ')
        for n in range(ndim):
            # Write binbounds per dim
            bb_file.write(f'{np.linspace(minlist[n], maxlist[n], nbins_per_dim[n] + 1)}\t')
        # Min/Max pcoord
        bb_file.write(f'\nLagging pcoord in each dimension: {minlist}\n')
        bb_file.write(f'Leading pcoord in each dimension: {maxlist}\n')
        # Bottlenecks bins exist
        if bottleneck:
            bb_file.write(f'Number of {strict_Z=} bottleneck bins filled: {n_bottleneck_filled} / {max_bottleneck}\n')
            for n in active_dims:
                if direction[n] in [0, 1, 86]:
                    bb_file.write(f'Dimension {n} forward bottleneck walker at: {[bottlenecks_forward[n]]}\n')
                if direction[n] in [0, -1, 86]:
                    bb_file.write(f'Dimension {n} backward bottleneck walker at: {[bottlenecks_reverse[n]]}\n')
            bb_file.write('\n')
        else:
            bb_file.write('\n')

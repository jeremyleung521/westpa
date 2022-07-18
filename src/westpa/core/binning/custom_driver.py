import logging
import math
import operator

import numpy as np

from westpa.core.we_driver import WEDriver

log = logging.getLogger(__name__)


def _sort_walkers_identity(we_driver, ibin, status, **kwargs):
    '''A function that, given  sorts the walkers based on a given criteria. Status indicate which method it's from. The int
    arguments mean the following:

    status = 0    _split_by_weight()
    status = 1    _merge_by_weight()
    status = 2    _adjust_count()
    status = 3    _sort_by_threshold()
    status = 4    _merge_by_threshold()
    '''
    log.debug('using we_driver._sort_walkers_identity')
    segments = np.array(sorted(ibin, key=operator.attrgetter('weight')), dtype=np.object_)
    weights = np.array(list(map(operator.attrgetter('weight'), segments)))
    cumul_weight = 0

    #    if status == 0:
    #        ordered_array = segments[weights > we_driver.weight_split_threshold * ideal_weight]
    #    elif status == 1:
    #        cumul_weight = np.add.accumulate(weights)
    #        ordered_array = segments[cumul_weight <= ideal_weight * self.weight_merge_cutoff]
    if status == 2:
        ordered_array = segments
    elif status == 3:
        ordered_array = segments[weights > we_driver.largest_allowed_weight]
    elif status == 4:
        ordered_array = segments[weights < we_driver.smallest_allowed_weight]
    else:
        print("I don\'t know what happened")

    print(status)
    return segments, weights, ordered_array, cumul_weight


class CustomDriver(WEDriver):
    #    def _split_by_weight(self, bin, target_count, ideal_weight):
    #        '''Split overweight particles'''
    #
    #        if len(bin) > 0:
    #            assert target_count > 0
    #        print(self.sorting_function_kwargs)
    #        kwargs = {ideal_weight=ideal_weight}
    #        segments, weights, to_split, _ = self.sorting_function(self, bin, 0, **kwargs)
    #
    #        for segment in to_split:
    #            m = int(math.ceil(segment.weight / ideal_weight))
    #            bin.remove(segment)
    #            new_segments_list = self._split_walker(segment, m, bin)
    #            bin.update(new_segments_list)
    #
    #    def _merge_by_weight(self, bin, target_count, ideal_weight):
    #        '''Merge underweight particles'''
    #
    #        while True:
    #            segments, weights, to_merge, cumul_weight = self.sorting_function(self, bin, 1)
    #
    #            if len(to_merge) < 2:
    #                return
    #            bin.difference_update(to_merge)
    #            new_segment, parent = self._merge_walkers(to_merge, cumul_weight, bin)
    #            bin.add(new_segment)
    #
    def _adjust_count(self, bin, subgroups, target_count):
        weight_getter = operator.attrgetter('weight')
        # Order subgroups by the sum of their weights.
        if len(subgroups) > target_count:
            sorted_subgroups = [set()]
            for i in bin:
                sorted_subgroups[0].add(i)
        else:
            sorted_subgroups = sorted(subgroups, key=lambda gp: sum(seg.weight for seg in gp))
        # Loops over the groups, splitting/merging until the proper count has been reached.  This way, no trajectories are accidentally destroyed.

        # split
        while len(bin) < target_count:
            for i in sorted_subgroups:
                log.debug('adjusting counts by splitting')
                # always split the highest probability walker into two
                segments, _, _, _ = self.sorting_function(self, bin, 2)
                bin.remove(segments[-1])
                i.remove(segments[-1])
                new_segments_list = self._split_walker(segments[-1], 2, bin)
                i.update(new_segments_list)
                bin.update(new_segments_list)

                if len(bin) == target_count:
                    break

        # merge
        while len(bin) > target_count:
            sorted_subgroups.reverse()
            # Adjust to go from lowest weight group to highest to merge
            for i in sorted_subgroups:
                # Ensures that there are least two walkers to merge
                if len(i) > 1:
                    log.debug('adjusting counts by merging')
                    # always merge the two lowest-probability walkers
                    segments, _, _, _ = self.sorting_function(self, bin, 2)
                    segments = sorted(i, key=weight_getter)
                    bin.difference_update(segments[:2])
                    i.difference_update(segments[:2])
                    merged_segment, parent = self._merge_walkers(segments[:2], cumul_weight=None, bin=bin)
                    i.add(merged_segment)
                    bin.add(merged_segment)

                    # As long as we're changing the merge_walkers and split_walkers, adjust them so that they don't update the bin within the function
                    # and instead update the bin here.  Assuming nothing else relies on those.  Make sure with grin.
                    # in bash, "find . -name \*.py | xargs fgrep -n '_merge_walkers'"
                    if len(bin) == target_count:
                        break

    def _split_by_threshold(self, bin, subgroup):
        # split to satisfy weight thresholds
        # this splits walkers that are too big
        segments, weights, to_split, _ = self.sorting_function(self, bin, 4)

        for segment in to_split:
            m = int(math.ceil(segment.weight / self.largest_allowed_weight))
            bin.remove(segment)
            subgroup.remove(segment)
            new_segments_list = self._split_walker(segment, m, bin)
            bin.update(new_segments_list)
            subgroup.update(new_segments_list)

    def _merge_by_threshold(self, bin, subgroup):
        # merge to satisfy weight thresholds
        # this gets rid of weights that are too small
        while True:
            segments, weights, to_merge, cumul_weight = self.sorting_function(self, bin, 5)

            if len(to_merge) < 2:
                return
            bin.difference_update(to_merge)
            subgroup.difference_update(to_merge)
            new_segment, parent = self._merge_walkers(to_merge, cumul_weight, bin)
            bin.add(new_segment)
            subgroup.add(new_segment)

    def __init__(self):
        self.sorting_function = _sort_walkers_identity
        self.sorting_function_kwargs = {}
        super().__init__()
        # super().__init__(rc=None, system=None)


#    def assign(self, segments, initializing=False):
#        '''Assign segments to initial and final bins, and update the (internal) lists of used and available
#        initial states. This function is adapted to the MAB scheme, so that the inital and final segments are
#        sent to the bin mapper at the same time, otherwise the inital and final bin boundaries can be inconsistent.'''
#
#        log.debug("CustomDriver in use.")
#        # collect initial and final coordinates into one place
#        n_segments = len(segments)
#        all_pcoords = np.empty((n_segments * 2, self.system.pcoord_ndim + 2), dtype=self.system.pcoord_dtype)
#
#        for iseg, segment in enumerate(segments):
#            all_pcoords[iseg] = np.append(segment.pcoord[0, :], [segment.weight, 0.0])
#            all_pcoords[n_segments + iseg] = np.append(segment.pcoord[-1, :], [segment.weight, 1.0])
#
#        # assign based on initial and final progress coordinates
#        assignments = self.bin_mapper.assign(all_pcoords)
#        initial_assignments = assignments[:n_segments]
#        if initializing:
#            final_assignments = initial_assignments
#        else:
#            final_assignments = assignments[n_segments:]
#
#        initial_binning = self.initial_binning
#        final_binning = self.final_binning
#        flux_matrix = self.flux_matrix
#        transition_matrix = self.transition_matrix
#        for (segment, iidx, fidx) in zip(segments, initial_assignments, final_assignments):
#            initial_binning[iidx].add(segment)
#            final_binning[fidx].add(segment)
#            flux_matrix[iidx, fidx] += segment.weight
#            transition_matrix[iidx, fidx] += 1
#
#        n_recycled_total = self.n_recycled_segs
#        n_new_states = n_recycled_total - len(self.avail_initial_states)
#
#        log.debug(
#            '{} walkers scheduled for recycling, {} initial states available'.format(
#                n_recycled_total, len(self.avail_initial_states)
#            )
#        )
#
#        if n_new_states > 0:
#            return n_new_states
#        else:
#            return 0

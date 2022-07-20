import logging
import math
import operator

import numpy as np

from westpa.core.we_driver import WEDriver

log = logging.getLogger(__name__)


def _sort_walkers_identity(we_driver, ibin, status, **kwargs):
    '''A function that, given  sorts the walkers based on a given criteria. Status indicate which method it's from. The int
    arguments mean the following:

    status = 0    _run_we() - not doing any sorting
    status = 1    _split_by_weight() - check upper ideal weight threshold
    status = 2    _merge_by_weight() - check lower ideal weight threshold
    status = 3    _adjust_count()
    status = 4    _split_by_threshold() - check upper weight threshold
    status = 5    _merge_by_threshold() - check lower weight threshold
    status = 6    _run_we() - merging all segs in one group
    '''
    log.debug('using we_driver._sort_walkers_identity')
    segments = np.array(sorted(ibin, key=operator.attrgetter('weight')), dtype=np.object_)
    weights = np.array(list(map(operator.attrgetter('weight'), segments)))
    cumul_weight = 0

    if status == 0:  # run_we - not doing any sorting
        ordered_array = []
    elif status == 1:  # _split_by_weight() - check upper ideal weight threshold
        ideal_weight = kwargs['ideal_weight']
        ordered_array = segments[weights > we_driver.weight_split_threshold * ideal_weight]
    elif status == 2:  # _merge_by_weight() - check lower ideal weight threshold
        cumul_weight = np.add.accumulate(weights)
        ideal_weight = kwargs['ideal_weight']
        ordered_array = segments[cumul_weight <= ideal_weight * we_driver.weight_merge_cutoff]
    elif status == 3:  # _adjust_count()
        ordered_array = segments
    elif status == 4:  # _split_by_threshold() - check upper weight threshold
        ordered_array = segments[
            weights > we_driver.largest_allowed_weight and not np.isclose(weights, we_driver.largest_allowed_weight)
        ]
    elif status == 5:  # _merge_by_threshold() - check lower weight threshold
        ordered_array = segments[
            weights < we_driver.smallest_allowed_weight and not np.isclose(weights, we_driver.smallest_allowed_weight)
        ]
    elif status == 6:  # _run_we - merging all segs in one group
        ordered_array = np.add.accumulate(weights)
    else:
        print(status)
        print("Not sure why this is triggered")

    return segments, weights, ordered_array, cumul_weight


class CustomDriver(WEDriver):
    def _split_by_weight(self, bin, target_count, ideal_weight):
        '''Split overweight particles'''
        if len(bin) > 0:
            assert target_count > 0
        self.sorting_function_kwargs['ideal_weight'] = ideal_weight
        segments, weights, to_split, _ = self.sorting_function(self, bin, 1, **self.sorting_function_kwargs)
        if self.sorting_function_kwargs['scheme'] == 'list':
            for segment in to_split:
                m = int(math.ceil(segment.weight / ideal_weight))
                bin.remove(segment)
                new_segments_list = self._split_walker(segment, m, bin)
                bin.update(new_segments_list)
        elif self.sorting_function_kwargs['scheme'] == 'paired':
            for segment in to_split:
                pass

    def _merge_by_weight(self, bin, target_count, ideal_weight):
        '''Merge underweight particles'''
        self.sorting_function_kwargs['ideal_weight'] = ideal_weight
        while True:
            segments, weights, to_merge, cumul_weight = self.sorting_function(self, bin, 2, **self.sorting_function_kwargs)
            if len(to_merge) < 2:
                return
            bin.difference_update(to_merge)
            new_segment, parent = self._merge_walkers(to_merge, cumul_weight, bin)
            bin.add(new_segment)

    def adjust_count(self, bin, subgroups, target_count):
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
                segments, _, _, _ = self.sorting_function(self, bin, 3, **self.sorting_function_kwargs)
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
                    segments, _, _, _ = self.sorting_function(self, bin, 3, **self.sorting_function_kwargs)
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
        segments, weights, to_split, _ = self.sorting_function(self, bin, 4, **self.sorting_function_kwargs)

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
            segments, weights, to_merge, cumul_weight = self.sorting_function(self, bin, 5, **self.sorting_function_kwargs)

            if len(to_merge) < 2:
                return
            bin.difference_update(to_merge)
            subgroup.difference_update(to_merge)
            new_segment, parent = self._merge_walkers(to_merge, cumul_weight, bin)
            bin.add(new_segment)
            subgroup.add(new_segment)

    def _run_we(self):
        '''Run recycle/split/merge. Do not call this function directly; instead, use
        populate_initial(), rebin_current(), or construct_next().'''
        self._recycle_walkers()

        # sanity check
        self._check_pre()

        # Regardless of current particle count, always split overweight particles and merge underweight particles
        # Then and only then adjust for correct particle count
        total_number_of_subgroups = 0
        total_number_of_particles = 0
        for (ibin, bin) in enumerate(self.next_iter_binning):
            if len(bin) == 0:
                continue

            # Splits the bin into subgroups as defined by the called function
            target_count = self.bin_target_counts[ibin]
            subgroups = self.subgroup_function(self, ibin, **self.subgroup_function_kwargs)
            total_number_of_subgroups += len(subgroups)
            # Clear the bin
            segments, weights, _, _ = self.sorting_function(self, bin, 0, **self.sorting_function_kwargs)
            ideal_weight = weights.sum() / target_count
            bin.clear()
            # Determines to see whether we have more sub bins than we have target walkers in a bin (or equal to), and then uses
            # different logic to deal with those cases.  Should devolve to the Huber/Kim algorithm in the case of few subgroups.
            if len(subgroups) >= target_count:
                for i in subgroups:
                    # Merges all members of set i.  Checks to see whether there are any to merge.
                    if len(i) > 1:
                        _, _, to_merge, _ = self.sorting_function(self, i, 6, **self.sorting_function_kwargs)
                        (segment, parent) = self._merge_walkers(
                            list(i),
                            to_merge,
                            i,
                        )
                        i.clear()
                        i.add(segment)
                    # Add all members of the set i to the bin.  This keeps the bins in sync for the adjustment step.
                    bin.update(i)

                if len(subgroups) > target_count:
                    self._adjust_count(bin, subgroups, target_count)

            if len(subgroups) < target_count:
                for i in subgroups:
                    self._split_by_weight(i, target_count, ideal_weight)
                    self._merge_by_weight(i, target_count, ideal_weight)
                    # Same logic here.
                    bin.update(i)
                if self.do_adjust_counts:
                    # A modified adjustment routine is necessary to ensure we don't unnecessarily destroy trajectory pathways.
                    self._adjust_count(bin, subgroups, target_count)
            if self.do_thresholds:
                for i in subgroups:
                    self._split_by_threshold(bin, i)
                    self._merge_by_threshold(bin, i)
                for iseg in bin:
                    if iseg.weight > self.largest_allowed_weight or iseg.weight < self.smallest_allowed_weight:
                        log.warning(
                            f'Unable to fulfill threshold conditions for {iseg}. The given threshold range is likely too small.'
                        )
            total_number_of_particles += len(bin)
        log.debug('Total number of subgroups: {!r}'.format(total_number_of_subgroups))

        self._check_post()

        self.new_weights = self.new_weights or []

        log.debug('used initial states: {!r}'.format(self.used_initial_states))
        log.debug('available initial states: {!r}'.format(self.avail_initial_states))

    def __init__(self, rc=None, system=None):
        self.sorting_function = _sort_walkers_identity
        self.sorting_function_kwargs = {}
        super().__init__()

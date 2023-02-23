from typing import List


class Range:
    start: int
    end: int

    def __init__(self, start, end):
        self.start = start
        self.end = end


def union_sets(sets: List[List[Range]]) -> List[Range]:
    """
    Returns the union of all sets in the input list, represented by a list of Range objects.
    Each set is represented by a list of Range objects.

    Parameters:
    sets: a list of sets of points, where each set is represented by a list of Exon objects,
    and each object represents a range.

    Returns:
    A list of Exon objects, where each object represents a range in the union of all sets.
    """
    result = []
    current_start = None
    current_end = None
    for r in sorted([r for set_ in sets for r in set_], key=lambda r: r.start):
            if current_end is None:
                current_start = r.start
                current_end = r.end
            elif r.start > current_end:
                result.append(Range(current_start, current_end))
                current_start = r.start
                current_end = r.end
            else:
                current_end = r.end
    return result


def divide_set(annotations: List[Range], distance: int) -> List[List[Range]]:
    """
  Divides a set of points into smaller sets such that no annotation in any set is more than
  `distance` away from all other annotations in that set. The input set is not modified.

  Parameters:
  annotations: a list of ExonAnnotation objects.
  distance: the maximum distance between annotations in a set.

  Returns:
  A list of lists of ExonAnnotation objects, where each list represents a smaller set, and each object
  represents an annotation in that set.
  """
    # sort the annotations in the set by their start values
    annotations = sorted(annotations, key=lambda a: a.start)

    # initialize an empty list to hold the smaller sets
    smaller_sets = []

    # initialize a last_end value to be the end value of the first annotation in the set
    last_end = annotations[0].end

    # initialize a list to hold the annotations in the current smaller set
    current_set = []

    # iterate over the annotations in the set
    for a in annotations:
        # if the annotation is more than `distance` away from the last_end value, add the current smaller set to the list of smaller sets
        # and start a new smaller set with the current annotation
        if a.start - last_end > distance:
            smaller_sets.append(current_set)
            current_set = [a]
            last_end = a.end
        # otherwise, add the annotation to the current smaller set
        else:
            current_set.append(a)
            last_end = a.end

    # add the final smaller set to the list of smaller sets
    smaller_sets.append(current_set)

    # return the list of smaller sets
    return smaller_sets
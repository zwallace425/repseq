# Imports
from typing import List, Dict, Set, Optional
from tqdm import tqdm


score_memo: Dict[str, int] = {}
def score(
    sequence1: str,
    sequence2: str,
    mr: int = 1,
    mp: int = 1,
    ip: int = 1,
) -> int:
    """
    Given two sequences, a match reward mr, a mismatch penalty mp, and an
    indel penalty ip, compute the alignment score between the two sequences.
    The alignment score between two sequences is computed as (
         number of matches × match reward
         - number of mismatches × mismatch penalty
         - number of indels × indel penalty
    )

    :param sequence1: <str> First sequence.
    :param sequence2: <str> Second sequence.
    :param mr: <int> Match reward. Default: 1.
    :param mp: <int> Mismatch penalty. Default: 1.
    :param ip: <int> Indel penalty. Default: 1.

    :return: <int> Score between the two sequences.

    Examples
    --------
    >>> score('GATACA-', '-ATAGAT')
    1
    >>> score('GATACA-', '-ATAGAT', 1, 1, 1)
    1
    >>> score('GATACA-', '-ATAGAT', 2, 1, 2)
    3
    >>> score('GATACA-', '-ATAGAT', 0, 1, 1)
    -3
    """
    # Check if input in cache
    input1 = f"{sequence1}_{sequence2}_{mr}_{mp}_{ip}"
    if input1 in score_memo:
        return score_memo[input1]
    input2 = f"{sequence2}_{sequence1}_{mr}_{mp}_{ip}"
    if input2 in score_memo:
        return score_memo[input2]

    # Compute score
    value = 0
    for a, b in zip(sequence1, sequence2):
        if a == "-" and b == "-":
            continue
        elif a == "-" or b == "-":
            value -= ip
        elif a == b:
            value += mr
        else:
            value -= mp

    # Store result in cache and return
    score_memo[input1] = value
    return value


global_score_memo: Dict[str, int] = {}
def global_alignment_score(
    sequence1: str,
    sequence2: str,
    mr: int = 1,
    mp: int = 1,
    ip: int = 1,
) -> int:
    """
    Given two sequences, a match reward mr, a mismatch penalty mp, and an
    indel penalty ip, perform global alignment between the two sequences
    and return the score.

    :param sequence1: <str> First sequence.
    :param sequence2: <str> Second sequence.
    :param mr: <int> Match reward. Default: 1.
    :param mp: <int> Mismatch penalty. Default: 1.
    :param ip: <int> Indel penalty. Default: 1.

    :return: <int> Score between the two sequences.

    Examples
    --------
    >>> global_alignment_score('GATACA', 'ATAGAT')
    1
    >>> global_alignment_score('GATACA', 'ATAGAT', 1, 1, 1)
    1
    >>> global_alignment_score('GATACA', 'ATAGAT', 2, 1, 2)
    3
    >>> global_alignment_score('GATACA', 'ATAGAT', 0, 1, 1)
    -3
    >>> global_alignment_score('GATACA', 'ATAGATA', 0, 1, 1)
    -4
    """
    # Memoization
    input1 = f"{sequence1}_{sequence2}_{mr}_{mp}_{ip}"
    if input1 in global_score_memo:
        return global_score_memo[input1]
    input2 = f"{sequence2}_{sequence1}_{mr}_{mp}_{ip}"
    if input2 in global_score_memo:
        return global_score_memo[input2]

    # Initialization
    values = {i: {0: -ip * i} for i in range(len(sequence1) + 1)}
    values[0] = {j: -ip * j for j in range(len(sequence2) + 1)}

    # Algorithm
    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            if sequence1[i - 1] == sequence2[j - 1]:
                match = mr
            else:
                match = -mp
            v = max(
                values[i - 1][j] - ip,
                values[i][j - 1] - ip,
                values[i - 1][j - 1] + match,
            )
            values[i][j] = v
    value = values[i][j]
    global_score_memo[input1] = value
    return value


def greedy_representative_subset(
    S: List[str],
    d: Optional[int] = None,
    mr: int = 1,
    mp: int = 1,
    ip: int = 1,
    aligned: bool = True,
) -> Set[str]:
    """
    Given a set S of strings, a parameter delta d, a match reward, a mismatch penalty,
    and an indel penalty, compute the smallest subset s of S such that for each sequence
    in S, there exists a sequence in s whose alignment score is d or higher using a
    suboptimal algorithm.

    Algorithm:
    1. Iterating over all the sequences.
        1.1 Add the sequence to the subset if there is no sequence in the subset with
            an alignment score greater or equal than the parameter delta.
        1.2 If there is such sequence, add that later sequence in a 'keep' list.

    2. Iterate over the elements in the subset that are not in the 'keep' list and
       check if it is safe to remove them.

    :param S: <List[str]> Input set of sequences S.
    :param d: <int> Delta. Minimum alignment score between sequences. If the alignemnt
                    score between to sequences is greater or equal than d we say that
                    the sequences represent each other. Default: 75% of the length of
                    the input sequences.
    :param mr: <int> Match reward. Default: 1.
    :param mp: <int> Mismatch penalty. Default: 1.
    :param ip: <int> Indel penalty. Default: 1.
    :param aligned: <bool> boolean Flag. Indicates whether if the input set S is aligned
                           or not. Default: True.

    :return: <List[str]> Subset s.
    """
    # Initialize parameter d
    if d is None:
        l_s = len(S[0])
        d = l_s * 75 // 100

    # Initialize variables
    s: Set[str] = set()
    keep: Set[str] = set()

    # Main algorithm
    for i, sequence1 in enumerate(tqdm(S)):
        for sequence2 in keep:
            if aligned:
                alignment_score = score(sequence1, sequence2, mr, mp, ip)
            else:
                alignment_score = global_alignment_score(sequence1, sequence2, mr, mp, ip)
            if alignment_score >= d:
                break
        else:
            for sequence2 in s - keep:
                if aligned:
                    alignment_score = score(sequence1, sequence2, mr, mp, ip)
                else:
                    alignment_score = global_alignment_score(sequence1, sequence2, mr, mp, ip)
                if alignment_score >= d:
                    keep.add(sequence1)
                    s.add(sequence1)
                    s.remove(sequence2)
                    break
            else:
                s.add(sequence1)

    # Prune unnecessary sequences
    for sequence1 in s - keep:
        for sequence2 in keep:
            if aligned:
                alignment_score = score(sequence1, sequence2, mr, mp, ip)
            else:
                alignment_score = global_alignment_score(sequence1, sequence2, mr, mp, ip)
            if alignment_score >= d:
                s.remove(sequence1)
                break

    return s


def greedy_representative_subset_v2(
    S: List[str],
    d: Optional[int] = None,
    mr: int = 1,
    mp: int = 1,
    ip: int = 1,
    aligned: bool = True,
) -> Set[str]:
    """
    Given a set S of strings, a parameter delta d, a match reward, a mismatch penalty,
    and an indel penalty, compute the smallest subset s of S such that for each sequence
    in S, there exists a sequence in s whose alignment score is d or higher using a
    suboptimal algorithm.

    Algorithm:
    1. Compute a n x n cover matrix where element i,j is equal to 1 if the global
       alignment score between sequences i and j is >= than d, 0 otherwise.
    2. While the subset cover is different to the identity vector:
        2.1 Add to subset the sequence that makes the sum of the subset cover highest.
    Note: The subset cover is computed as the elementwise 'or' operation between all
    elements of the subset.

    :param S: <List[str]> Input set of sequences S.
    :param d: <int> Delta. Minimum alignment score between sequences. If the alignemnt
                    score between to sequences is greater or equal than d we say that
                    the sequences represent each other. Default: 75% of the length of
                    the input sequences.
    :param mr: <int> Match reward. Default: 1.
    :param mp: <int> Mismatch penalty. Default: 1.
    :param ip: <int> Indel penalty. Default: 1.
    :param aligned: <bool> boolean Flag. Indicates whether if the input set S is aligned
                           or not. Default: True.

    :return: <Set[str]> Subset s.
    """
    # Initialize parameter d
    if d is None:
        l_s = len(S[0])
        d = l_s * 75 // 100

    # Compute al paired scores
    cover_matrix: List[List[int]] = []
    for i, sequence1 in enumerate(tqdm(S)):
        row = []
        for j in range(i):
            row.append(cover_matrix[j][i])
        row.append(1)
        for sequence2 in S[i + 1:]:
            if aligned:
                alignment_score = score(sequence1, sequence2, mr, mp, ip)
            else:
                alignment_score = global_alignment_score(sequence1, sequence2, mr, mp, ip)
            covered = 1 if alignment_score >= d else 0
            row.append(covered)
        cover_matrix.append(row)

    # Main algorithm
    subset_cover = [0] * len(cover_matrix[0])
    s: Set[str] = set()
    max_cover = 0
    while not all(subset_cover):
        element = None
        for i, element_cover in enumerate(cover_matrix):
            if S[i] in s:
                continue
            cover = sum(a | b for a, b in zip(subset_cover, element_cover))
            if cover > max_cover:
                max_cover = cover
                element = i
        subset_cover = [a | b for a, b in zip(subset_cover, cover_matrix[element])]
        s.add(S[element])

    return s

# Imports
from typing import List


score_memo = {}
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
    # TODO: add tests
    # TODO: add logging outputs
    # Check if input in cache
    input = f"{sequence1}_{sequence2}_{mr}_{mp}_{ip}"
    if input in score_memo:
        return score_memo[input]
    input2 = f"{sequence1}_{sequence2}_{mr}_{mp}_{ip}"
    if input2 in score_memo:
        return score_memo[input]

    # Compute score
    score = 0
    for a, b in zip(sequence1, sequence2):
        if a == "-" and b == "-":
            continue
        elif a == "-" or b == "-":
            score -= ip
        elif a == b:
            score += mr
        else:
            score -= mp

    # Store result in cache and return
    score_memo[input] = score
    return score


def representative_subset(
    S: List[str],
    d: int = None,
    mr: int = 1,
    mp: int = 1,
    ip: int = 1,
) -> List[str]:
    """
    Given a set S of strings, a parameter delta d, a match reward, a mismatch penalty, and an indel penalty,
    compute the smallest subset s of S such that for each sequence in S, there exists a sequence in s whose alignment score is d or higher.

    :param S: <List[str]> Input set of sequences S.
    :param d: <int> Delta. Minimum alignment score between sequences. If the alignemnt
                    score between to sequences is greater or equal than d we say that
                    the sequences represent each other. Default: 75% of the length of the input sequences.
    :param mr: <int> Match reward. Default: 1.
    :param mp: <int> Mismatch penalty. Default: 1.
    :param ip: <int> Indel penalty. Default: 1.

    :return: <List[str]> Subset s.
    """
    # TODO: add tests
    # Initialize parameter d
    if d == None:
        l_s = len(S[0])
        d = l_s * 75 // 100

    # Initialize variables
    s = set()
    keep = set()

    # Main algorithm
    for sequence in S:
        for sequence2 in keep:
            alignment_score = score(sequence, sequence2, mr, mp, ip)
            if alignment_score >= d:
                break
        else:
            for sequence2 in s-keep:
                alignment_score = score(sequence, sequence2, mr, mp, ip)
                if alignment_score >= d:
                    keep.add(sequence2)
                    break
            else:
                s.add(sequence)

    # Prune unnecessary sequences
    for sequence in s-keep:
        for sequence2 in keep:
            alignment_score = score(sequence, sequence2, mr, mp, ip)
            if alignment_score >= d:
                s.remove(sequence)
                break

    return s


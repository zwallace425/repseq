## repseq: Selecting Represenative Sequences from a Large Set

The purpose of the project is to devise an algorithm that will select that smallest
subset of Influenza A sequences from a large set of sequences for a given subtype 
that best represents the diversity within entire subtype.

### Comutational Problem

We are designing an algorithm based on the following computational problem statement:

`Definitions`:

Given a globally aligned set of genome strings S and a threshold parameter delta d, we wish to find the smallest subset s from S in which for any string within S, there exists a string in s such that their alignment score is at least d.

We can assume that any pair of sequences of the set S are globally aligned. Given a match reward, a mismatch penalty, and an indel penalty value, we can compute the alignment score between two aligned sequences as (number of matches * match reward - number of mismatches * mismatch penalty - number of indels * indel penalty).

`Input`: A set S of strings that have been aligned globally, a parameter delta d, a match reward, a mismatch penalty, and an indel penalty.

`Output`: The smallest subset s of S such that for each sequence in S, there exists a sequence in s whose alignment score is d or higher.

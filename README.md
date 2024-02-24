## repseq: Selecting Represenative Sequences from a Large Set

The purpose of the project is to devise an algorithm that will select that smallest
subset of Influenza A sequences from a large set of sequences for a given subtype 
that best represents the diversity within entire subtype.

### Comutational Problem

We are designing an algorithm based on the following computational problem statement:

Given strings v and w and scoring parameters, we define score(v,w) as the score of the global alignment between these strings. Given an algined string-set S, we define score score(S,w) as the maximal score(v,w) among all strings v in S.  

We can assume that any pair of sequences of the set S are globally aligned. Given a match reward, a mismatch penalty, and an indel penalty value, we can compute the global alignment score between two aligned sequences as (number of matches * match reward - number of mismatches * mismatch penalty - number of indels * indel penalty).

`Input`: A string-set set S, a parameter delta d, and alignment scoring parameters (match reward, mismatch penalty, and indel penalty). 

`Output`: The smallest subset S’ of S such that for each sequence w in S, score(S’,w)>=d.

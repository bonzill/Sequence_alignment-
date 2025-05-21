def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    # Initialize the scoring matrix
    n = len(seq1) + 1
    m = len(seq2) + 1
    score_matrix = [[0] * m for _ in range(n)]
    
    # Initialize the first row and column
    for i in range(n):
        score_matrix[i][0] = gap_penalty * i
    for j in range(m):
        score_matrix[0][j] = gap_penalty * j

    # Fill in the scoring matrix
    for i in range(1, n):
        for j in range(1, m):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback to find the optimal alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n - 1, m - 1

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        if i > 0 and j > 0 and current_score == score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and current_score == score_matrix[i - 1][j] + gap_penalty:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), score_matrix[-1][-1]

# Sample sequences (DNA sequences)
seq1 = "AGCTGACGTA"
seq2 = "AGCTGATGCA"

# Note: These sequences represent segments of DNA that may be compared for similarity
# in genetic studies, which can help identify evolutionary relationships or functional 
# similarities between genes.

# Perform alignment
aligned_seq1, aligned_seq2, score = needleman_wunsch(seq1, seq2)

print("Aligned Sequences:")
print(aligned_seq1)
print(aligned_seq2)
print("Alignment Score:", score)

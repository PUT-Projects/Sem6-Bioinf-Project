def calculate_hash(s, base, powers):
    hash_value = 0
    for i in range(len(s)):
        hash_value += (ord(s[i]) - ord('A') + 1) * powers[len(s) - i - 1]
    return hash_value

def find_offset(s1, s2):
    n = len(s1)
    if n != len(s2):
        return -1  # Ensure that both strings are of equal length

    base = 5
    powers = [1] * n

    # Compute powers of base up to n-1
    for i in range(1, n):
        powers[i] = powers[i - 1] * base

    # Hash for the first string
    hash_s1 = calculate_hash(s1[1:], base, powers)
    current_hash_s2 = calculate_hash(s2[:-1], base, powers)  # Initial hash for first 'n' chars of s2

    if current_hash_s2 == hash_s1:
        return 1  # Offset found

    for i in range(2, n - 1):
        # Update hash for s1 (remove first char)
        hash_s1 -= (ord(s1[i - 2]) - ord('A') + 1) * powers[n - i - 1]

        # Update hash for s2 (remove char from end)
        current_hash_s2 -= (ord(s2[i - 2]) - ord('A') + 1) * powers[n - 2]




        if hash_s1 == current_hash_s2:
            return i

    return -1  # No offset found


def find_offset_naive(s1, s2):
    n = len(s1)
    if n != len(s2):
        return -1  # Ensure that both strings are of equal length

    for i in range(1, n):
        if s1[i:] == s2[:-i]:
            return i

    return -1  # No offset found

print(find_offset_naive("AACTGC", "ACTGCG"))  # Should return 1
print(find_offset_naive("AACTGC", "CTGCCT"))  # Should return 3
print(find_offset_naive("AACTGC", "TGCCTA"))  # Should return 3
print(find_offset("AACTGC", "TGCCTA"))  # Should return 3
print(find_offset("AACTGC", "GAATTC"))  # Should return 6 (full length as no overlap)
print(find_offset("AACTGC", "AACTGC"))  # Should return 0
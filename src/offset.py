def find_offset_naive(s1: str, s2: str) -> int:

    n = len(s1)
    if n != len(s2):
        return -1  # Ensure that both strings are of equal length

    for i in range(1, n):
        if s1[i:] == s2[:-i]:
            return i

    return n  # No offset found
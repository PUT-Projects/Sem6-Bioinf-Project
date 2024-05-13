def find_offset_direct(str1, str2):
    len1 = len(str1)
    len2 = len(str2)

    # Start from index 1 because the first characters are different
    for i in range(1, len1 - 1):
        # Check if str2 fits within the remaining length of str1 from index i
        if str1[i:] == str2[:len1 - i]:
            return i

    return len1  # If no match is found, return the length of str1

# Examples
print(find_offset_direct("AACTGC", "ACTCCG"))  # Example 1: Offset should be 1
print(find_offset_direct("AACTGC", "TGCCTA"))  # Example 2: Offset should be 3
print(find_offset_direct("AACTGC", "GAATTC"))  # Example 3: No overlap; should return 6
print(find_offset_direct("AACTGC", "AACTGC"))  # Example 4: Offset should be 0 (even if they differ at index 0 in typical inputs)

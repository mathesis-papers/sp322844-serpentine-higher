import hashlib
import itertools
import json
import math
from fractions import Fraction


def matrix(size):
    return [
        [
            column + size * (row - 1) if row % 2 else size * row - column + 1
            for column in range(1, size + 1)
        ]
        for row in range(1, size + 1)
    ]


def determinant(entries):
    size = len(entries)
    if size == 0:
        return 1
    work = [list(row) for row in entries]
    previous = 1
    sign = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (row for row in range(pivot_index, size) if work[row][pivot_index]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    pivot * work[row][column]
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                assert numerator % previous == 0
                work[row][column] = numerator // previous
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def rank(entries):
    if not entries:
        return 0
    work = [[Fraction(value) for value in row] for row in entries]
    found = 0
    for column in range(len(work[0])):
        pivot_row = next(
            (row for row in range(found, len(work)) if work[row][column]), None
        )
        if pivot_row is None:
            continue
        work[found], work[pivot_row] = work[pivot_row], work[found]
        pivot = work[found][column]
        work[found] = [value / pivot for value in work[found]]
        for row in range(found + 1, len(work)):
            multiplier = work[row][column]
            work[row] = [
                value - multiplier * normalized
                for value, normalized in zip(work[row], work[found])
            ]
        found += 1
        if found == len(work):
            break
    return found


def multiply(left, right):
    columns = list(zip(*right))
    return [
        [sum(first * second for first, second in zip(row, column)) for column in columns]
        for row in left
    ]


def characteristic_coefficients(entries):
    size = len(entries)
    power = [[int(row == column) for column in range(size)] for row in range(size)]
    traces = []
    for exponent in range(1, size + 1):
        power = multiply(power, entries)
        traces.append(sum(power[index][index] for index in range(size)))
    coefficients = [Fraction(1)]
    for degree in range(1, size + 1):
        coefficients.append(
            -sum(coefficients[degree - index] * traces[index - 1] for index in range(1, degree + 1))
            / degree
        )
    assert all(value.denominator == 1 for value in coefficients)
    return [int(value) for value in coefficients]


def quantities(size):
    parity = size % 2
    trace = size * (size * size + parity) // 2
    magnitude = Fraction(size * size * (3 * (1 + size * size) - 2 * (2 + size * size) * parity), 12)
    assert magnitude.denominator == 1
    return trace, int(magnitude)


def polynomial_product(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for first_degree, first_value in enumerate(left):
        for second_degree, second_value in enumerate(right):
            result[first_degree + second_degree] += first_value * second_value
    return result


def exponential_coefficient(polynomial, degree, sign=1):
    return sum(
        Fraction(value * sign ** (degree - index), math.factorial(degree - index))
        for index, value in enumerate(polynomial) if index <= degree
    )


def permanent(entries):
    size = len(entries)
    states = [0] * (1 << size)
    states[0] = 1
    for subset in range(1 << size):
        row = subset.bit_count()
        if row == size:
            continue
        for column in range(size):
            if not subset & (1 << column):
                states[subset | (1 << column)] += states[subset] * entries[row][column]
    return states[-1]


def aggregate_sums(size):
    parity = size % 2
    upper = Fraction(size * size * (size + 1) * (2 * size + 1) + 3 * size * (size + parity), 12)
    lower = Fraction(size * size * (size * size + 1), 2) + quantities(size)[0] - upper
    return upper, lower


def linked_quantities(size):
    return Fraction(size * size * (size * size - 1), 12), Fraction(size * size * (size * size + 1), 2)


def check_siblings():
    permanent_values = []
    for size in range(10):
        row_coefficients = [Fraction(1)]
        column_coefficients = [Fraction(1)]
        for index in range(1, size + 1):
            shifted = Fraction(2 * size * index - size + 1, 2)
            sign = (-1) ** (index + 1)
            centered = Fraction(2 * index - size - 1, 2)
            row_coefficients = polynomial_product(row_coefficients, [shifted, sign])
            column_coefficients = polynomial_product(column_coefficients, [1, centered])
        assert all(not column_coefficients[index] for index in range(1, size + 1, 2))
        formula = sum(math.factorial(index) * math.factorial(size - index)
                      * row_coefficients[index] * column_coefficients[index]
                      for index in range(size + 1))
        actual = permanent(matrix(size))
        assert formula == actual
        permanent_values.append(actual)
    assert permanent_values[1:7] == [1, 11, 490, 60916, 15745548, 7477647372]
    for size in range(2, 17):
        pyramidal, magic_total = linked_quantities(size)
        assert pyramidal.denominator == magic_total.denominator == 1
        for sign in [1, -1]:
            entries = [[row + sign * column for column in range(1, size + 1)] for row in range(1, size + 1)]
            expected_trace = size * (size + 1) if sign == 1 else 0
            assert characteristic_coefficients(entries) == [1, -expected_trace, -sign * pyramidal] + [0] * (size - 2)
            assert rank(entries) == 2
        assert sum((row - column) ** 2 for row in range(1, size + 1) for column in range(1, size + 1)) == 2 * pyramidal
    for size in range(2, 17):
        entries = matrix(size)
        trace, magnitude = quantities(size)
        previous = entries
        current = multiply(entries, entries)
        for exponent in range(3, 7):
            following = multiply(current, entries)
            assert following == [[trace * current[row][column] + magnitude * previous[row][column]
                                  for column in range(size)] for row in range(size)]
            previous, current = current, following
    for size in range(1, 33):
        entries = matrix(size)
        upper = sum(entries[row][column] for row in range(size) for column in range(row, size))
        lower = sum(entries[row][column] for row in range(size) for column in range(row + 1))
        assert (upper, lower) == aggregate_sums(size)
        assert entries[-1][0] == (size * size if size % 2 == 0 else size * size - size + 1)
    denominator = polynomial_product([(-1) ** index * math.comb(5, index) for index in range(6)], [1, 2, 1])
    upper_numerator = [0, 1, 3, 7, 3, 2]
    lower_numerator = [0, 1, 5, 13, 9, 4]
    even_numerator = [0, 5, 43, 43, 5]
    odd_numerator = [0, 6, 20, 6]
    for degree in range(65):
        magnitude = quantities(degree)[1]
        corrected = (exponential_coefficient([0, 3, 15, 12, 2], degree)
                     + exponential_coefficient([0, -3, 9, -6, 1], degree, -1)) / 12
        assert corrected == Fraction(magnitude, math.factorial(degree))
        shifted = (exponential_coefficient([24, 21, 9, 2, 1], degree + 2, -1)
                   + exponential_coefficient([36, -33, 15, -4, 2], degree + 2)) / 12
        assert shifted == Fraction(magnitude, math.factorial(degree + 2))
        even_value = quantities(2 * degree)[1]
        odd_value = quantities(2 * degree + 1)[1]
        assert exponential_coefficient([0, 5, 29, 24, 4], degree) == Fraction(even_value, math.factorial(degree))
        assert exponential_coefficient([0, 6, 19, Fraction(32, 3), Fraction(4, 3)], degree) == Fraction(odd_value, math.factorial(degree))
        factorial_inverse = lambda index: Fraction(1, math.factorial(index)) if index >= 0 else Fraction(0)
        even_weighted = (8 * factorial_inverse(2 * degree - 2) + factorial_inverse(2 * degree - 4)
                         + 2 * factorial_inverse(2 * degree - 1) + 6 * factorial_inverse(2 * degree - 3)) / 4
        odd_weighted = (6 * factorial_inverse(2 * degree - 2) + 6 * factorial_inverse(2 * degree - 1)
                        + factorial_inverse(2 * degree - 3)) / 12
        assert even_weighted == Fraction(even_value, math.factorial(2 * degree))
        assert odd_weighted == Fraction(odd_value, math.factorial(2 * degree + 1))
        for parity, numerator in [(0, even_numerator), (1, odd_numerator)]:
            coefficient = sum((-1) ** index * math.comb(5, index) * quantities(2 * (degree - index) + parity)[1]
                              for index in range(6) if index <= degree)
            assert coefficient == (numerator[degree] if degree < len(numerator) else 0)
        for aggregate_index, numerator in [(0, upper_numerator), (1, lower_numerator)]:
            coefficient = sum(denominator[index] * aggregate_sums(degree - index)[aggregate_index]
                              for index in range(len(denominator)) if index <= degree)
            assert coefficient == (numerator[degree] if degree < len(numerator) else 0)
        lower_egf = (exponential_coefficient([0, 3], degree, -1)
                     + exponential_coefficient([0, 21, 78, 54, 8], degree)) / 24
        assert lower_egf == aggregate_sums(degree)[1] / math.factorial(degree)
        pyramidal = linked_quantities(degree)[0]
        assert exponential_coefficient([0, 0, Fraction(1, 2), Fraction(1, 2), Fraction(1, 12)], degree) == pyramidal / math.factorial(degree)
        pyramidal_numerator = [0, 0, 1, 1]
        coefficient = sum((-1) ** index * math.comb(5, index) * linked_quantities(degree - index)[0]
                          for index in range(6) if index <= degree)
        assert coefficient == (pyramidal_numerator[degree] if degree < len(pyramidal_numerator) else 0)
    assert Fraction(5, 24) != Fraction(5, 2)
    assert 3 + 5 + 13 * 2 == 34 != aggregate_sums(2)[1]
    literal_lower_egf = (exponential_coefficient([0, 3], 2, -1)
                         + exponential_coefficient([0, 21, 102, 54], 2)) / 24
    assert literal_lower_egf * math.factorial(2) == 10 != aggregate_sums(2)[1]
    for index in range(1, 65):
        assert Fraction(1, quantities(2 * index)[1]) == Fraction(1, index * index) - Fraction(4, 4 * index * index + 1)
        assert Fraction(1, quantities(2 * index + 1)[1]) == Fraction(3, index * (index + 1)) - Fraction(12, (2 * index + 1) ** 2)
        pyramidal, magic_total = linked_quantities(index)
        assert 1 / magic_total == 2 * (Fraction(1, index * index) - Fraction(1, index * index + 1))
        if index >= 2:
            assert 1 / pyramidal == 12 * (Fraction(1, index * index - 1) - Fraction(1, index * index))
    terms = 4096
    even_partial = math.fsum(1 / quantities(2 * index)[1] for index in range(1, terms + 1))
    odd_partial = math.fsum(1 / quantities(2 * index + 1)[1] for index in range(1, terms + 1))
    even_constant = 2 + math.pi ** 2 / 6 - math.pi / math.tanh(math.pi / 2)
    odd_constant = 15 - 3 * math.pi ** 2 / 2
    assert abs(even_constant - even_partial) <= 1 / (12 * terms ** 3) + 1e-14
    assert abs(odd_constant - odd_partial) <= 1 / (4 * terms ** 3) + 1e-14
    pyramidal_partial = math.fsum(1 / float(linked_quantities(index)[0]) for index in range(2, terms + 1))
    pyramidal_alternating = math.fsum((-1) ** index / float(linked_quantities(index)[0]) for index in range(2, terms + 1))
    magic_partial = math.fsum(1 / float(linked_quantities(index)[1]) for index in range(1, terms + 1))
    magic_alternating = math.fsum((-1) ** (index + 1) / float(linked_quantities(index)[1]) for index in range(1, terms + 1))
    assert abs(pyramidal_partial - (21 - 2 * math.pi ** 2)) <= 16 / (3 * terms ** 3) + 1e-14
    assert abs(pyramidal_alternating - (math.pi ** 2 - 9)) <= 16 / (3 * terms ** 3) + 1e-14
    assert abs(magic_partial - (1 + math.pi ** 2 / 3 - math.pi / math.tanh(math.pi))) <= 2 / (3 * terms ** 3) + 1e-14
    assert abs(magic_alternating - (math.pi / math.sinh(math.pi) + math.pi ** 2 / 6 - 1)) <= 2 / (3 * terms ** 3) + 1e-14
    return {
        'permanent': {'sizes': [0, 9], 'values': permanent_values, 'method': 'subset dynamic programming compared with symmetric-polynomial formula'},
        'ordinary_power_recurrence': {'sizes': [2, 16], 'powers': [3, 6]},
        'a002415_matrix_families': {'sizes': [2, 16], 'entries': ['i+j', 'i-j'], 'characteristic_polynomials_and_rank': 'pass', 'squared_difference_formula': 'pass'},
        'a002415_a037270_reciprocal_series': {'partial_fraction_indices_through': 64, 'terms_each': terms, 'series': 4, 'absolute_roundoff_allowance': 1e-14, 'pyramidal_absolute_tail_bound': 16 / (3 * terms ** 3), 'magic_total_absolute_tail_bound': 2 / (3 * terms ** 3)},
        'diagonal_aggregates_and_last_entry': {'sizes': [1, 32]},
        'scalar_and_aggregate_generating_functions': {'coefficients_through_degree': 64, 'arithmetic': 'exact rational coefficients'},
        'reciprocal_partial_fractions': {'indices': [1, 64], 'arithmetic': 'exact rational identities'},
        'reciprocal_constant_numerical_checks': {'terms_each': terms, 'absolute_tolerance': 1e-14, 'even_tail_upper_bound': 1 / (12 * terms ** 3), 'odd_tail_upper_bound': 1 / (4 * terms ** 3)},
        'conventional_egf_correction': {'degree_two_source_shifted_weight': '5/24', 'degree_two_standard_weight': '5/2'},
        'a325656_literal_products': {'size': 2, 'ordinary_literal': 34, 'exponential_literal': 10, 'correct_value': 8},
    }


def main():
    if not __debug__:
        raise RuntimeError('Run without Python optimization so that every exact check is enabled.')
    minor_count = 0
    minor_digest = hashlib.sha256()
    for size in range(11):
        entries = matrix(size)
        trace, magnitude = quantities(size)
        for degree in range(size + 1):
            total = 0
            for subset in itertools.combinations(range(size), degree):
                value = determinant([[entries[row][column] for column in subset] for row in subset])
                if degree >= 3:
                    assert value == 0
                total += value
                minor_count += 1
                minor_digest.update(f'{size}:{subset}:{value}\n'.encode())
            expected = 1 if degree == 0 else trace if degree == 1 else -magnitude if degree == 2 else 0
            assert total == expected, (size, degree, total, expected)
    characteristic_digest = hashlib.sha256()
    for size in range(1, 17):
        actual = characteristic_coefficients(matrix(size))
        trace, magnitude = quantities(size)
        expected = [1, -trace] if size == 1 else [1, -trace, -magnitude] + [0] * (size - 2)
        assert actual == expected, (size, actual, expected)
        characteristic_digest.update(f'{size}:{actual}\n'.encode())
    diagonal_count = 0
    for size in range(33):
        entries = matrix(size)
        trace, magnitude = quantities(size)
        assert rank(entries) == min(size, 2)
        assert determinant(entries) == (1 if size <= 1 else -5 if size == 2 else 0)
        centered = [Fraction(2 * index - size - 1, 2) for index in range(1, size + 1)]
        shifted = [Fraction(2 * size * index - size + 1, 2) for index in range(1, size + 1)]
        signs = [(-1) ** (index + 1) for index in range(1, size + 1)]
        for row in range(size):
            assert sum(entries[row]) == size * size * (row + 1) - size * (size - 1) // 2
            for column in range(size):
                assert entries[row][column] == shifted[row] + signs[row] * centered[column]
            source_row = row + 1
            antidiagonal = entries[row][size - source_row]
            assert antidiagonal == ((size - 1) * source_row + 1 if source_row % 2 else (size + 1) * source_row - size)
        for column in range(size):
            expected = Fraction(size * (size * size + 1), 2) + size % 2 * centered[column]
            assert sum(entries[row][column] for row in range(size)) == expected
        for offset in range(size):
            length = size - offset
            parity = length % 2
            upper = sum(entries[index][index + offset] for index in range(length))
            lower = sum(entries[index + offset][index] for index in range(length))
            assert upper == Fraction(size * (2 * length * length + 1 - (-1) ** length), 4)
            assert lower == Fraction(
                size * length * (size + offset) + length * (1 - (-1) ** offset)
                + (-1) ** offset * (size - 2 * offset) * parity, 2
            )
            source_lower = Fraction(
                2 * (-1 + (-1) ** size) * offset - 2 * offset * offset * size
                + size * (2 - (-1) ** offset - (-1) ** size + 2 * size * size), 4
            )
            assert lower == source_lower
            diagonal_count += 2
        if size >= 2:
            assert magnitude > 0
    assert [quantities(size)[1] for size in range(1, 7)] == [0, 5, 6, 68, 50, 333]
    assert [matrix(4)[index][3 - index] for index in range(4)] == [4, 6, 10, 16]
    assert [matrix(4)[index][3 - index] for index in range(4)] != [6, 8, 10, 12]
    assert {index * (2 * index * index - 2 * index + 1) % 5 for index in range(5)} == {0, 1, 4}
    assert {index * index * (2 * index + 1) % 5 for index in range(5)} == {0, 3, 4}
    expected_numerator = [0, 0, 5, 6, 43, 20, 43, 6, 5, 0]
    for degree in range(65):
        coefficient = sum(
            (-1) ** index * math.comb(5, index) * quantities(degree - 2 * index)[1]
            for index in range(6) if degree >= 2 * index
        )
        assert coefficient == (expected_numerator[degree] if degree < 10 else 0)
    result = {
        'arithmetic': 'exact integers and fractions',
        'principal_minors': {'sizes': [0, 10], 'count': minor_count, 'sha256': minor_digest.hexdigest()},
        'characteristic_polynomials': {'sizes': [1, 16], 'sha256': characteristic_digest.hexdigest()},
        'factorization_rank_row_column_sums': {'sizes': [0, 32]},
        'diagonal_sums': {'sizes': [1, 32], 'count': diagonal_count},
        'antidiagonal_entries': {'sizes': [1, 32], 'size_four': [4, 6, 10, 16], 'a074147_row_four': [6, 8, 10, 12]},
        'ordinary_generating_function': {'coefficients_through_degree': 64, 'numerator_coefficients': expected_numerator},
        'superdiagonal_halved_mod_five': {'even_size': [0, 1, 4], 'odd_size': [0, 3, 4]},
        'coverage': 'Finite checks of the explicit matrices and polynomial identities; the all-size proofs are in note.tex.',
        'sibling_checks': check_siblings(),
        'status': 'pass',
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

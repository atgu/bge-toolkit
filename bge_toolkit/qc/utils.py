from enum import Enum

import hail as hl

ALL_GROUP = 'ALL'
ALL_AGG = 'ALL'


class JoinType(str, Enum):
    """
    The join type to use.

    Attributes:
        outer: Outer join.
        inner: Inner join.
    """

    outer = 'outer'
    inner = 'inner'


def full_outer_join(left: hl.MatrixTable, right: hl.MatrixTable, join_type: JoinType) -> hl.MatrixTable:
    left_sample_counter = left.aggregate_cols(hl.agg.counter(left.col_key[0]))
    right_sample_counter = right.aggregate_cols(hl.agg.counter(right.col_key[0]))

    left_bad = [f'{k!r}: {v}' for k, v in left_sample_counter.items() if v > 1]
    right_bad = [f'{k!r}: {v}' for k, v in right_sample_counter.items() if v > 1]
    if left_bad or right_bad:
        raise ValueError(
            f"Found duplicate sample IDs:\n" f"  left:  {', '.join(left_bad)}\n" f"  right: {', '.join(right_bad)}"
        )

    included = set(left_sample_counter.keys()).intersection(set(right_sample_counter.keys()))

    lit = hl.literal(included, dtype=hl.tset(hl.tstr))
    left = left.filter_cols(lit.contains(left.col_key[0]))
    right = right.filter_cols(lit.contains(right.col_key[0]))

    joined = hl.experimental.full_outer_join_mt(left, right)

    if join_type == JoinType.inner:
        return joined.filter_rows(hl.is_defined(joined.left_row) & hl.is_defined(joined.right_row))

    return joined

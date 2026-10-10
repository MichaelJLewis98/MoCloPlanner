from types import SimpleNamespace

import pytest

from moclo_assembly import get_valid_parts_for_position
from moclo_standards import (
    CIDAR,
    get_expected_flanks,
    validate_part_chain,
)


def vector(code="renamed-vector", family="DVK", start="A", end="E"):
    return SimpleNamespace(
        code=code,
        backbone=family,
        five_prime_overhang=start,
        three_prime_overhang=end,
    )


def part(code, five_prime, three_prime, part_type="part"):
    return SimpleNamespace(
        code=code,
        five_prime_overhang=five_prime,
        three_prime_overhang=three_prime,
        part_type=part_type,
        ori="",
        contains=code,
    )


def test_cidar_l1_path_uses_explicit_vector_flanks_not_vector_name():
    destination = vector(code="DVK_AE")

    assert [get_expected_flanks(destination, index, "L1") for index in range(4)] == [
        ("A", "B"),
        ("B", "C"),
        ("C", "D"),
        ("D", "E"),
    ]
    assert get_expected_flanks(vector(code="a-renamed-vector"), 0, "L1") == ("A", "B")


def test_cidar_l1_internal_flanks_are_fixed_when_vector_boundaries_change():
    destination = vector(code="DVK_HI", start="H", end="I")

    assert [get_expected_flanks(destination, index, "L1") for index in range(4)] == [
        ("H", "B"),
        ("B", "C"),
        ("C", "D"),
        ("D", "I"),
    ]


def test_cidar_l2_uses_fixed_efg_internal_flanks_between_vector_boundaries():
    destination = vector(code="arbitrary-name", family="DVC", start="A", end="H")

    assert [get_expected_flanks(destination, index, "L2") for index in range(4)] == [
        ("A", "E"),
        ("E", "F"),
        ("F", "G"),
        ("G", "H"),
    ]
    assert get_expected_flanks(destination, 4, "L2") == (None, None)


def test_vector_boundary_cannot_repeat_an_internal_flank():
    destination = vector(code="DVK_BE", start="B", end="E")

    with pytest.raises(ValueError, match="reuse a flank"):
        get_expected_flanks(destination, 0, "L1")


def test_vector_without_explicit_flanks_is_rejected_even_if_name_has_range():
    destination = vector(code="DVK_AE", start="", end="")

    with pytest.raises(ValueError, match="explicit single-letter"):
        get_expected_flanks(destination, 0, "L1")


def test_part_chain_accepts_ordered_flanks_and_rejects_incompatible_part():
    destination = vector()
    valid_chain = [
        part("promoter", "A", "B"),
        part("rbs", "B", "C"),
        part("cds", "C", "D"),
        part("terminator", "D", "E"),
    ]
    validate_part_chain(destination, valid_chain, "L1", CIDAR)

    incompatible_chain = valid_chain[:2] + [part("wrong-cds", "B", "D")] + valid_chain[3:]
    with pytest.raises(ValueError, match="position 3 requires C->D"):
        validate_part_chain(destination, incompatible_chain, "L1", CIDAR)


def test_part_type_does_not_override_incompatible_flanks():
    destination = vector()
    wrong_promoter = part("promoter", "C", "D", part_type="Promoter")

    _, _, candidates = get_valid_parts_for_position(
        [wrong_promoter], destination, 0, "L1"
    )

    assert candidates == []


def test_incomplete_part_chain_fails_validation():
    with pytest.raises(ValueError, match="path is incomplete"):
        validate_part_chain(vector(), [part("promoter", "A", "B")], "L1")


def test_l0_manual_parts_chain_explicit_junctions_to_vector_boundaries():
    destination = vector(code="label-only", family="DVA", start="F", end="G")
    inserts = [part("fragment-1", "F", "Q"), part("fragment-2", "Q", "G")]

    validate_part_chain(destination, inserts, "L0")

    with pytest.raises(ValueError, match="preceding flank Q"):
        validate_part_chain(
            destination,
            [part("fragment-1", "F", "Q"), part("fragment-2", "R", "G")],
            "L0",
        )

"""Part filtering and labeling helpers independent of the Streamlit interface."""

from moclo_standards import DEFAULT_STANDARD, get_expected_flanks


def clean(value):
    return str(value).strip().upper()


def get_next_insert_index(selected_parts):
    for idx, part in enumerate(selected_parts):
        if part is None:
            return idx
    return len(selected_parts)


def has_matching_overhangs(part, expected_5, expected_3):
    five_prime = clean(getattr(part, "five_prime_overhang", ""))
    three_prime = clean(getattr(part, "three_prime_overhang", ""))
    return (five_prime, three_prime) == (expected_5, expected_3)


def get_valid_parts_for_position(
    parts,
    backbone,
    position_index,
    assembly_level,
    show_all_parts=False,
    standard=DEFAULT_STANDARD,
):
    expected_5, expected_3 = get_expected_flanks(
        backbone, position_index, assembly_level, standard
    )
    if expected_5 is None or expected_3 is None:
        return expected_5, expected_3, []

    if show_all_parts:
        return expected_5, expected_3, list(parts)

    valid_parts = [
        part for part in parts
        if has_matching_overhangs(part, expected_5, expected_3)
    ]
    return expected_5, expected_3, valid_parts


def part_identifier(part):
    fields = [
        getattr(part, "code", ""),
        getattr(part, "inventory_sheet", ""),
        getattr(part, "location", ""),
        getattr(part, "ori", ""),
        getattr(part, "contains", ""),
        getattr(part, "part_type", ""),
        getattr(part, "five_prime_overhang", ""),
        getattr(part, "three_prime_overhang", ""),
    ]
    return " | ".join(str(field).strip() for field in fields)


def dedupe_parts(parts):
    unique_parts = []
    seen = set()
    for part in parts:
        token = part_identifier(part)
        if token in seen:
            continue
        seen.add(token)
        unique_parts.append(part)
    return unique_parts


def chain_signature(backbone, selected_parts):
    backbone_code = clean(part_identifier(backbone) or "backbone")
    part_codes = [
        clean(part_identifier(part) or f"slot_{index}")
        for index, part in enumerate(selected_parts)
        if part is not None
    ]
    return "__".join([backbone_code] + part_codes)


def construct_label(backbone, inserts):
    backbone_code = str(getattr(backbone, "code", "")).strip()
    insert_codes = [str(getattr(insert, "code", "")).strip() for insert in inserts]
    return " / ".join([backbone_code] + insert_codes)

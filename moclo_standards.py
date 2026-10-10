"""Assembly syntax and vector constraints for supported MoClo standards."""

from dataclasses import dataclass
from string import ascii_uppercase
from typing import Mapping


@dataclass(frozen=True)
class AssemblyStandard:
    """Overhang-path rules for an assembly standard, independent of reaction mix."""

    name: str
    vector_families_by_level: Mapping[str, tuple[str, ...]]
    inventory_part_levels: Mapping[str, str]
    manual_insert_levels: tuple[str, ...]
    variant_assembly_levels: tuple[str, ...]
    level1_internal_flanks: tuple[str, ...]
    level2_internal_flanks: tuple[str, ...]

    def allowed_vector_families(self, assembly_level: str) -> set[str]:
        return set(self.vector_families_by_level.get(assembly_level.upper(), ()))

    def inventory_part_level(self, assembly_level: str) -> str:
        return self.inventory_part_levels.get(assembly_level.upper(), "")

    def uses_manual_inserts(self, assembly_level: str) -> bool:
        return assembly_level.upper() in self.manual_insert_levels

    def supports_variant_assemblies(self, assembly_level: str) -> bool:
        return assembly_level.upper() in self.variant_assembly_levels

    @property
    def assembly_levels(self) -> tuple[str, ...]:
        return tuple(self.vector_families_by_level)


CIDAR = AssemblyStandard(
    name="CIDAR",
    vector_families_by_level={
        "L0": ("DVA",),
        "L1": ("DVK",),
        "L2": ("DVC", "DVA", "DVS"),
    },
    inventory_part_levels={"L0": "L0", "L1": "L0", "L2": "L1"},
    manual_insert_levels=("L0",),
    variant_assembly_levels=("L1", "L2"),
    level1_internal_flanks=("B", "C", "D"),
    level2_internal_flanks=("E", "F", "G"),
)

STANDARDS: dict[str, AssemblyStandard] = {CIDAR.name: CIDAR}
DEFAULT_STANDARD = CIDAR
allowed_backbones = {
    level: DEFAULT_STANDARD.allowed_vector_families(level)
    for level in DEFAULT_STANDARD.vector_families_by_level
}


def _vector_boundaries(vector) -> tuple[str, str]:
    code = getattr(vector, "code", "<unknown>")
    start = str(getattr(vector, "five_prime_overhang", "") or "").strip().upper()
    end = str(getattr(vector, "three_prime_overhang", "") or "").strip().upper()
    if len(start) != 1 or len(end) != 1 or start not in ascii_uppercase or end not in ascii_uppercase:
        raise ValueError(
            f"Vector {code} must have explicit single-letter 5' and 3' flank values in the inventory."
        )
    if start == end:
        raise ValueError(f"Vector {code} must have different 5' and 3' flanks.")
    return start, end


def _path_with_internal_flanks(
    vector,
    start: str,
    end: str,
    internal_flanks: tuple[str, ...],
    assembly_level: str,
    standard: AssemblyStandard,
) -> tuple[str, ...]:
    code = getattr(vector, "code", "<unknown>")
    path = (start, *internal_flanks, end)
    if any(flank not in ascii_uppercase for flank in internal_flanks):
        raise ValueError(
            f"{standard.name} {assembly_level} internal flanks must be single letters A-Z."
        )
    if len(set(path)) != len(path):
        joined_path = "->".join(path)
        raise ValueError(
            f"Vector {code} boundaries {start}->{end} reuse a flank in the "
            f"{standard.name} {assembly_level} path {joined_path}."
        )
    return path


def validate_vector(vector, assembly_level: str, standard: AssemblyStandard = DEFAULT_STANDARD) -> None:
    """Validate vector family and explicit boundary flanks against the standard."""
    code = getattr(vector, "code", "<unknown>")
    family = str(getattr(vector, "backbone", "") or "").strip().upper()
    allowed = standard.allowed_vector_families(assembly_level)
    if family not in allowed:
        choices = ", ".join(sorted(allowed)) or "none"
        raise ValueError(
            f"Vector {code} has backbone family {family or '<blank>'}; "
            f"{standard.name} {assembly_level} expects one of: {choices}."
        )
    start, end = _vector_boundaries(vector)
    if assembly_level.upper() == "L1":
        _path_with_internal_flanks(
            vector,
            start,
            end,
            standard.level1_internal_flanks,
            "L1",
            standard,
        )
    elif assembly_level.upper() == "L2":
        _path_with_internal_flanks(
            vector,
            start,
            end,
            standard.level2_internal_flanks,
            "L2",
            standard,
        )


def expected_flank_path(
    vector,
    assembly_level: str,
    standard: AssemblyStandard = DEFAULT_STANDARD,
) -> list[tuple[str, str]]:
    """Return the ordered insert flank pairs allowed by a vector and standard."""
    validate_vector(vector, assembly_level, standard)
    start, end = _vector_boundaries(vector)
    level = assembly_level.upper()

    if level == "L1":
        path = (start, *standard.level1_internal_flanks, end)
        return list(zip(path, path[1:]))

    if level == "L2":
        path = (start, *standard.level2_internal_flanks, end)
        return list(zip(path, path[1:]))

    if level in standard.manual_insert_levels:
        # L0 inputs can be a single long fragment or a variable chain. Validate
        # their explicit junctions directly rather than inventing intermediates.
        return []

    raise ValueError(f"{standard.name} does not define assembly rules for level {level}.")


def get_expected_flanks(
    vector,
    insert_index: int,
    assembly_level: str,
    standard: AssemblyStandard = DEFAULT_STANDARD,
) -> tuple[str | None, str | None]:
    """Get an insert's required 5'/3' flanks by its zero-based chain position."""
    if insert_index < 0:
        raise ValueError("Insert index must be zero or greater.")
    path = expected_flank_path(vector, assembly_level, standard)
    if insert_index >= len(path):
        return None, None
    return path[insert_index]


def validate_part_chain(
    vector,
    inserts,
    assembly_level: str,
    standard: AssemblyStandard = DEFAULT_STANDARD,
    *,
    require_complete: bool = True,
) -> None:
    """Reject missing, incompatible, or incomplete overhang chains."""
    path = expected_flank_path(vector, assembly_level, standard)
    if not inserts:
        raise ValueError("Add at least one insert before calculating a reaction.")

    if assembly_level.upper() in standard.manual_insert_levels:
        start, end = _vector_boundaries(vector)
        next_flank = start
        for index, part in enumerate(inserts):
            code = getattr(part, "code", getattr(part, "label", f"insert {index + 1}"))
            actual_5 = str(getattr(part, "five_prime_overhang", "") or "").strip().upper()
            actual_3 = str(getattr(part, "three_prime_overhang", "") or "").strip().upper()
            if not actual_5 or not actual_3:
                raise ValueError(f"Part {code} needs explicit 5' and 3' flank values.")
            if actual_5 != next_flank:
                raise ValueError(
                    f"Part {code} starts at {actual_5}; position {index + 1} "
                    f"must start at the preceding flank {next_flank}."
                )
            next_flank = actual_3
        if require_complete and next_flank != end:
            raise ValueError(
                f"The assembly chain ends at {next_flank}; vector {getattr(vector, 'code', '<unknown>')} "
                f"requires the chain to end at {end}."
            )
        return

    for index, part in enumerate(inserts):
        if index >= len(path):
            raise ValueError(
                f"The assembly has more inserts than the vector path supports "
                f"({len(path)} for {standard.name} {assembly_level})."
            )
        expected_5, expected_3 = path[index]
        code = getattr(part, "code", getattr(part, "label", f"insert {index + 1}"))
        actual_5 = str(getattr(part, "five_prime_overhang", "") or "").strip().upper()
        actual_3 = str(getattr(part, "three_prime_overhang", "") or "").strip().upper()
        if not actual_5 or not actual_3:
            raise ValueError(f"Part {code} needs explicit 5' and 3' flank values.")
        if (actual_5, actual_3) != (expected_5, expected_3):
            raise ValueError(
                f"Part {code} has flanks {actual_5}->{actual_3}; position {index + 1} "
                f"requires {expected_5}->{expected_3} for {standard.name} {assembly_level}."
            )

    if require_complete and len(inserts) != len(path):
        next_5, next_3 = path[len(inserts)]
        raise ValueError(
            f"The assembly path is incomplete: {len(inserts)} of {len(path)} parts were added; "
            f"the next required flanks are {next_5}->{next_3}."
        )

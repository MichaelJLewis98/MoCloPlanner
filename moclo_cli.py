"""Interactive command-line interface for MoCloPlanner."""

import argparse
from pathlib import Path
from typing import List

from moclo_core import (
    DNAItem, Inventory, ReactionProtocol, build_level0_reaction_rows,
    build_reaction, calculate_mq_water, create_manual_dna_item,
    reaction_to_dataframe,
)
from moclo_standards import (
    DEFAULT_STANDARD,
    STANDARDS,
    validate_part_chain,
    validate_vector,
)

def choose_inventory_match(matches: List[DNAItem]) -> DNAItem | None:
    if len(matches) == 1:
        return matches[0]
    if not matches:
        return None

    print("Multiple stock rows match this code:")
    for index, match in enumerate(matches, start=1):
        print(f"  {index}. {match.label}")
    while True:
        choice = input("Choose a stock row number, or press Enter to re-enter the code: ").strip()
        if not choice:
            return None
        try:
            index = int(choice)
        except ValueError:
            index = 0
        if 1 <= index <= len(matches):
            return matches[index - 1]
        print(f"Enter a number from 1 to {len(matches)}.")


def choose_item_by_code(inventory: Inventory, prompt: str) -> DNAItem:
    while True:
        code = input(prompt).strip()
        item = choose_inventory_match(inventory.find_all_by_code(code))
        if item is not None:
            return item

        if inventory.find_all_by_code(code):
            continue
        print(f"No exact match for '{code}'.")
        matches = inventory.search(code)
        if matches:
            print("Possible matches:")
            for match in matches[:10]:
                print(f" - {match.label}")
        else:
            print("No partial matches found.")

def prompt_for_manual_insert(default_level: str = "L0") -> DNAItem:
    print("\nAdding a manual DNA part. Enter values or press Enter to cancel.")
    code = input("Code (required): ").strip()
    if code == "":
        raise KeyboardInterrupt("Manual entry cancelled")

    contains = input("Contains / Name: ").strip() or code
    part_type = input("Part type [e.g. CDS, promoter, backbone]: ").strip() or ""

    # parse numeric fields with simple retries
    while True:
        length_raw = input("Length (bp): ").strip()
        try:
            length_bp = float(length_raw)
            break
        except Exception:
            print("Please enter a numeric length in base pairs.")

    while True:
        conc_raw = input("Concentration (ng/uL): ").strip()
        try:
            concentration_ng_ul = float(conc_raw)
            if concentration_ng_ul <= 0:
                raise ValueError
            break
        except Exception:
            print("Please enter a positive concentration (ng/uL).")

    five_prime_overhang = input("5' flank: ").strip().upper()
    three_prime_overhang = input("3' flank: ").strip().upper()

    level = default_level
    location = "manual"
    backbone = ""
    resistance = ""
    notes = "(manual entry)"

    return create_manual_dna_item(
        code=code,
        contains=contains,
        part_type=part_type,
        length_bp=length_bp,
        concentration_ng_ul=concentration_ng_ul,
        level=level,
        location=location,
        backbone=backbone,
        resistance=resistance,
        notes=notes,
        five_prime_overhang=five_prime_overhang,
        three_prime_overhang=three_prime_overhang,
    )

def choose_multiple_items(inventory: Inventory, default_level: str = "L0") -> List[DNAItem]:
    print(f"\nEnter insert codes one by one for {default_level}.")
    print("Press Enter on an empty line when done.")
    print("Type MANUAL to add a PCR fragment or other DNA part not present in the inventory.")

    selected: List[DNAItem] = []

    while True:
        code = input("Insert code: ").strip()
        if code == "":
            break

        if code.lower() == "manual":
            manual_item = prompt_for_manual_insert(default_level=default_level)
            selected.append(manual_item)
            print(f"Added manual part: {manual_item.label}")
            continue

        matches = inventory.find_all_by_code(code)
        if matches:
            matches = [
                match for match in matches
                if match.level.strip().upper() in {default_level.upper(), ""}
            ]
            item = choose_inventory_match(matches)
            if item is None:
                continue
        else:
            print(f"No exact match for '{code}'. Trying partial search...")
            matches = [match for match in inventory.search(code) if match.level == default_level or match.level == ""]
            if matches:
                for match in matches[:10]:
                    print(f" - {match.label}")
            else:
                print("No matches found.")
            continue

        if item.level and item.level.strip().upper() != default_level.upper():
            print(f"{item.code} is level {item.level}, but {default_level} parts are required here.")
            continue

        selected.append(item)
        print(f"Added: {item.label}")

    return selected

def choose_mode() -> str:
    print("\nChoose calculator mode:")
    print("  1 - General MoClo reaction plan")
    print("  2 - Level 0 insert/vector calculator template")

    while True:
        choice = input("Mode [1/2]: ").strip() or "1"
        if choice == "1":
            return "general"
        if choice == "2":
            return "level0"
        print("Please enter 1 or 2.")

def prompt_reaction_settings() -> ReactionProtocol:
    print("\nEnter fixed reaction settings.")
    total_reaction_ul = float(input("Total reaction volume in uL [7.5]: ").strip() or "7.5")
    buffer_ul = float(input("Buffer volume in uL [1.0]: ").strip() or "1.0")
    ligase_ul = float(input("Ligase volume in uL [0.5]: ").strip() or "0.5")
    enzyme_name = input("Restriction enzyme name [Restriction enzyme]: ").strip() or "Restriction enzyme"
    enzyme_ul = float(input("Restriction enzyme volume in uL [0.5]: ").strip() or "0.5")
    backbone_fmol = float(input("Backbone target in fmol [20]: ").strip() or "20")
    insert_ratio = float(input("Insert ratio [2]: ").strip() or "2")
    return ReactionProtocol(
        total_reaction_ul=total_reaction_ul,
        buffer_ul=buffer_ul,
        ligase_ul=ligase_ul,
        enzyme_ul=enzyme_ul,
        enzyme_name=enzyme_name,
        backbone_fmol=backbone_fmol,
        insert_ratio=insert_ratio,
    )

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Golden Gate / MoClo reaction planner")
    parser.add_argument("inventory_excel", type=Path, help="Path to the inventory Excel workbook")
    parser.add_argument("--sheet", default=None, help="Optional Excel sheet name")
    parser.add_argument(
        "--standard",
        choices=tuple(STANDARDS),
        default=DEFAULT_STANDARD.name,
        help="Assembly standard to use",
    )
    parser.add_argument("--backbone-fmol", type=float, default=20.0, help="Backbone target amount in fmol")
    parser.add_argument(
        "--insert-ratio",
        type=float,
        default=2.0,
        help="Per-insert ratio relative to backbone fmol",
    )
    parser.add_argument(
        "--total-reaction-ul",
        type=float,
        default=7.5,
        help="Total reaction volume in uL",
    )
    parser.add_argument(
        "--buffer-ul",
        type=float,
        default=1.0,
        help="Buffer volume in uL",
    )
    parser.add_argument(
        "--ligase-ul",
        type=float,
        default=0.5,
        help="Ligase volume in uL",
    )
    parser.add_argument(
        "--enzyme-ul",
        type=float,
        default=0.5,
        help="Restriction enzyme volume in uL",
    )
    parser.add_argument(
        "--enzyme-name",
        default="Restriction enzyme",
        help="Restriction enzyme name for the reaction protocol summary",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output Excel path for the reaction plan",
    )
    return parser.parse_args()

def choose_assembly_level(standard=DEFAULT_STANDARD) -> str:
    print("\nChoose assembly level:")
    levels = standard.assembly_levels
    for index, level in enumerate(levels):
        print(f"  {index} - {level}")

    while True:
        choice = input(f"Level [0-{len(levels) - 1}]: ").strip() or "0"
        if choice.isdigit() and int(choice) < len(levels):
            return levels[int(choice)]
        print(f"Please enter a number from 0 to {len(levels) - 1}.")


def choose_inventory_sheet_names(path: Path, assembly_level: str) -> List[str]:
    available = Inventory.available_inventory_sheet_names(path)
    if path.suffix.lower() == ".csv":
        return available

    default_sheet_map = {
        "L0": ["MoClo Dest Vectors", "MoClo Level 0 Parts"],
        "L1": ["MoClo Dest Vectors", "MoClo Level 0 Parts", "MoClo Level 1 TU's"],
        "L2": ["MoClo Dest Vectors", "MoClo Level 1 TU's", "MoClo Level 2 Multi TU's"],
    }

    suggested = [name for name in default_sheet_map.get(assembly_level, []) if name in available]

    print("\nAvailable inventory sheets:")
    for name in available:
        marker = "*" if name in suggested else " "
        print(f" {marker} {name}")

    if suggested:
        print("\nPress Enter to use the suggested sheets for this level.")
        print("Or type a comma-separated list of sheet names to override.")
        default_text = ", ".join(suggested)
        raw = input(f"Sheets [{default_text}]: ").strip()
        if raw == "":
            return suggested
    else:
        print("\nType a comma-separated list of sheet names to load.")
        raw = input("Sheets: ").strip()
        if raw == "":
            return available

    selected = [name.strip() for name in raw.split(",") if name.strip()]
    invalid = [name for name in selected if name not in available]
    if invalid:
        raise ValueError(
            "Unknown sheet name(s): " + ", ".join(invalid) + "\nAvailable sheets: " + ", ".join(available)
        )
    return selected

def main() -> int:
    args = parse_args()

    assembly_standard = STANDARDS[args.standard]
    assembly_level = choose_assembly_level(assembly_standard)
    selected_sheet_names = choose_inventory_sheet_names(args.inventory_excel, assembly_level)

    if args.inventory_excel.suffix.lower() == ".csv":
        inventory = Inventory.from_excel(args.inventory_excel, sheet_name=args.sheet)
    else:
        inventory = Inventory.from_excel_sheets(args.inventory_excel, selected_sheet_names)

    print(
        f"Loaded {len(inventory.items)} DNA items from {args.inventory_excel} "
        f"using sheet(s): {', '.join(selected_sheet_names)}"
    )
    print("\nExample inventory entries:")
    inventory.list_codes(limit=10)

    settings = ReactionProtocol(
        total_reaction_ul=args.total_reaction_ul,
        buffer_ul=args.buffer_ul,
        ligase_ul=args.ligase_ul,
        enzyme_ul=args.enzyme_ul,
        enzyme_name=args.enzyme_name,
        backbone_fmol=args.backbone_fmol,
        insert_ratio=args.insert_ratio,
    )

    mode = choose_mode()

    print(f"\nChoose destination/backbone vector for {assembly_level}")
    backbone = choose_item_by_code(inventory, "Backbone code: ")
    try:
        validate_vector(backbone, assembly_level, assembly_standard)
    except ValueError as exc:
        print(f"Invalid destination vector: {exc}")
        return 2

    print(f"\nChoose inserts or parts for {assembly_level}")
    required_part_level = assembly_standard.inventory_part_level(assembly_level)
    inserts = choose_multiple_items(inventory, default_level=required_part_level)
    if not inserts:
        print("No inserts selected. Exiting.")
        return 1
    try:
        validate_part_chain(backbone, inserts, assembly_level, assembly_standard)
    except ValueError as exc:
        print(f"Invalid assembly chain: {exc}")
        return 2

    if mode == "level0":
        level0_df = build_level0_reaction_rows(
            vector=backbone,
            inserts=inserts,
            fmol_target=settings.backbone_fmol,
            settings=settings,
        )

        print("\n=== Level 0 calculator table ===")
        print(level0_df.to_string(index=False))

        if args.output is not None:
            try:
                if args.output.suffix.lower() == ".csv":
                    level0_df.to_csv(args.output, index=False)
                else:
                    level0_df.to_excel(args.output, index=False)
                print(f"Saved reaction plan to {args.output}")
            except ImportError as exc:
                raise RuntimeError(
                    "Writing .xlsx files requires the 'openpyxl' package. "
                    "Install it with: pip install openpyxl\n"
                    "Or use an output filename ending in .csv instead."
                ) from exc
        return 0

    components = build_reaction(
        backbone=backbone,
        inserts=inserts,
        backbone_fmol=settings.backbone_fmol,
        insert_ratio=settings.insert_ratio,
    )
    reaction_df = reaction_to_dataframe(components)
    mq_water_ul = calculate_mq_water(components, settings)

    print("\n=== Reaction plan ===")
    print(reaction_df.to_string(index=False))
    total_dna_volume = reaction_df["Volume to add (uL)"].sum()
    print(f"\nTotal DNA volume: {total_dna_volume:.3f} uL")
    print(f"Fixed reagent volume: {settings.fixed_volume_ul:.3f} uL")
    print(f"Restriction enzyme: {settings.enzyme_name} ({settings.enzyme_ul:.3f} uL)")
    print(f"Buffer: {settings.buffer_ul:.3f} uL; ligase: {settings.ligase_ul:.3f} uL")
    print(f"MQ water: {mq_water_ul:.3f} uL")

    if args.output is not None:
        try:
            if args.output.suffix.lower() == ".csv":
                reaction_df.to_csv(args.output, index=False)
            else:
                reaction_df.to_excel(args.output, index=False)
            print(f"Saved reaction plan to {args.output}")
        except ImportError as exc:
            raise RuntimeError(
                "Writing .xlsx files requires the 'openpyxl' package. "
                "Install it with: pip install openpyxl\n"
                "Or use an output filename ending in .csv instead."
            ) from exc

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

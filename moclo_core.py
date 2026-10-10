"""Shared data models, inventory loading, and reaction calculations for MoCloPlanner."""

from dataclasses import dataclass
import math
from pathlib import Path
from typing import List, Optional

import pandas as pd

@dataclass
class DNAItem:
    code: str
    location: str
    level: str
    backbone: str
    resistance: str
    contains: str
    part_type: str
    length_bp: float
    concentration_ng_ul: float
    ori: str = ""
    a260_230: str = ""
    a260_280: str = ""
    sequenced: str = ""
    five_prime_overhang: str = ""
    three_prime_overhang: str = ""
    notes: str = ""
    inventory_sheet: str = ""

    @property
    def plasmid_backbone(self) -> str:
        return self.backbone

    @property
    def label(self) -> str:
        ori_part = f" | {self.ori}" if getattr(self, "ori", "") else ""
        try:
            length = float(self.length_bp)
            length_part = f"{length:.0f} bp" if math.isfinite(length) and length > 0 else "length missing"
        except (TypeError, ValueError):
            length_part = "length missing"
        stock_location = " / ".join(
            value for value in (self.inventory_sheet.strip(), self.location.strip()) if value
        )
        location_part = f" | {stock_location}" if stock_location else ""
        return f"{self.code}{ori_part} | {self.contains} | {self.part_type} | {length_part}{location_part}"


@dataclass
class CustomDNAInput:
    code: str
    contains: str
    part_type: str
    length_bp: float
    concentration_ng_ul: float
    level: str = "custom"
    location: str = "manual"
    backbone: str = ""
    resistance: str = ""
    notes: str = ""
    five_prime_overhang: str = ""
    three_prime_overhang: str = ""

    def to_dna_item(self) -> DNAItem:
        return DNAItem(
            code=self.code,
            location=self.location,
            level=self.level,
            backbone=self.backbone,
            resistance=self.resistance,
            contains=self.contains,
            part_type=self.part_type,
            length_bp=self.length_bp,
            concentration_ng_ul=self.concentration_ng_ul,
            five_prime_overhang=self.five_prime_overhang,
            three_prime_overhang=self.three_prime_overhang,
            notes=self.notes,
        )


@dataclass
class ReactionComponent:
    item: DNAItem
    target_fmol: float

    @property
    def required_ng(self) -> float:
        if not math.isfinite(self.target_fmol) or self.target_fmol <= 0:
            raise ValueError("Target DNA amount must be a positive finite value in fmol.")
        if not math.isfinite(self.item.length_bp) or self.item.length_bp <= 0:
            raise ValueError(f"Invalid DNA length for {_item_reference(self.item)}; expected positive bp.")
        # dsDNA approximation: 650 g/mol per bp, or 0.00065 ng per fmol-bp.
        return self.target_fmol * self.item.length_bp * 0.00065

    @property
    def required_volume_ul(self) -> float:
        if (
            not math.isfinite(self.item.concentration_ng_ul)
            or self.item.concentration_ng_ul <= 0
        ):
            raise ValueError(
                f"Invalid concentration for {_item_reference(self.item)}; "
                "expected a positive finite value in ng/uL."
            )
        return self.required_ng / self.item.concentration_ng_ul

    @property
    def role(self) -> str:
        if self.item.part_type.lower() in {"destination vector", "vector", "backbone"}:
            return "vector"
        return "insert"


@dataclass
class ReactionProtocol:
    total_reaction_ul: float = 7.5
    buffer_ul: float = 1.0
    ligase_ul: float = 0.5
    enzyme_ul: float = 0.5
    backbone_fmol: float = 20.0
    insert_ratio: float = 2.0
    enzyme_name: str = "Restriction enzyme"

    def __post_init__(self) -> None:
        positive = {
            "total reaction volume (uL)": self.total_reaction_ul,
            "backbone target (fmol)": self.backbone_fmol,
            "insert ratio": self.insert_ratio,
        }
        nonnegative = {
            "buffer volume (uL)": self.buffer_ul,
            "ligase volume (uL)": self.ligase_ul,
            "restriction enzyme volume (uL)": self.enzyme_ul,
        }
        for name, value in positive.items():
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name.capitalize()} must be a positive finite value.")
        for name, value in nonnegative.items():
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name.capitalize()} must be a nonnegative finite value.")
        if not isinstance(self.enzyme_name, str) or not self.enzyme_name.strip():
            raise ValueError("Restriction enzyme name cannot be blank.")

    @property
    def fixed_volume_ul(self) -> float:
        return self.buffer_ul + self.ligase_ul + self.enzyme_ul


# Backwards-compatible name for clients using the original class name.
ReactionSettings = ReactionProtocol


file_columns = {
    "Location" : "location",
    "Code" : "code",
    "Level" : "level",
    "Plasmid backbone" : "backbone",
    "Ori" : "ori",
    "Copy nr. (p/cell)" : "copy_number",
    "Resistance" : "resistance",
    "Contains (part/construct)" : "contains",
    "Part type" : "part_type",
    "5' flank" : "five_prime_overhang",
    "3' flank" : "three_prime_overhang",
    "Glyc. stock date" : "glyc_stock_date",
    "Description" : "description",
    "Depositor" : "depositor",
    "Plasmid conc. (ng/uL)" : "plasmid_concentration",
    "A260/230" : "a260_230",
    "A260/280" : "a260_280",
    "GC Content (%)" : "gc_content",
    "Length (bp)" : "length_bp",
    "Sequenced" : "sequenced",
    "Benchling file" : "benchling_file",
    "Notes": "notes"
}

# Average molecular weight per base pair (g / mol / bp).
# 650 is the standard approximation for double-stranded DNA.
BP_MW = 650.0

def _item_reference(item: DNAItem) -> str:
    stock_location = " / ".join(
        str(value).strip()
        for value in (
            getattr(item, "inventory_sheet", ""),
            getattr(item, "location", ""),
        )
        if value
    )
    code = str(getattr(item, "code", "<unknown>"))
    return f"{code} [{stock_location}]" if stock_location else code


class Inventory:
    @classmethod
    def from_excel_sheets(cls, path: str | Path, sheet_names: List[str]) -> "Inventory":
        path = Path(path)

        if path.suffix.lower() == ".csv":
            raise ValueError("Selecting multiple sheets is only supported for Excel workbooks, not CSV files.")

        all_frames: List[pd.DataFrame] = []
        for sheet_name in sheet_names:
            df = cls._read_excel_sheet_with_detected_header(path, sheet_name)
            df["inventory_sheet"] = sheet_name
            df = cls._normalise_dataframe(df)
            all_frames.append(df)

        if not all_frames:
            raise ValueError("No valid inventory sheets were loaded.")

        combined_df = pd.concat(all_frames, ignore_index=True)
        items = cls._rows_to_items(combined_df)
        return cls(items)
    def __init__(self, items: List[DNAItem]):
        seen: dict[tuple[str, str, str], DNAItem] = {}
        duplicates = []
        for item in items:
            key = tuple(
                value.strip().casefold()
                for value in (item.inventory_sheet, item.location, item.code)
            )
            if key in seen:
                duplicates.append(f"{item.code} [{item.inventory_sheet} / {item.location}]")
            else:
                seen[key] = item
        if duplicates:
            duplicate_text = ", ".join(sorted(set(duplicates), key=str.casefold))
            raise ValueError(
                "Duplicate inventory code/location record(s), ignoring case: "
                f"{duplicate_text}."
            )
        self.items = items

    @classmethod
    def from_excel(cls, path: str | Path, sheet_name: Optional[str] = None) -> "Inventory":
        path = Path(path)

        try:
            if path.suffix.lower() == ".csv":
                df = pd.read_csv(path)
                df = cls._normalise_dataframe(df)
                items = cls._rows_to_items(df)
                return cls(items)

            if sheet_name is not None:
                df = cls._read_excel_sheet_with_detected_header(path, sheet_name)
                df["inventory_sheet"] = sheet_name
                df = cls._normalise_dataframe(df)
                items = cls._rows_to_items(df)
                return cls(items)

            excel_file = pd.ExcelFile(path, engine="openpyxl")
            candidate_sheets = [
                name for name in excel_file.sheet_names
                if name.lower() not in {"general info", "sheet4"}
            ]

            all_frames: List[pd.DataFrame] = []
            for name in candidate_sheets:
                try:
                    sheet_df = cls._read_excel_sheet_with_detected_header(path, name)
                except ValueError:
                    continue
                sheet_df["inventory_sheet"] = name
                all_frames.append(cls._normalise_dataframe(sheet_df))

            if not all_frames:
                raise ValueError(
                    "Could not find a sheet with the expected inventory columns. "
                    "Try passing --sheet with one of the inventory sheet names."
                )

            combined_df = pd.concat(all_frames, ignore_index=True)
            items = cls._rows_to_items(combined_df)
            return cls(items)
        except ImportError as exc:
            raise RuntimeError(
                "Reading .xlsx files requires the 'openpyxl' package. "
                "Install it with: pip install openpyxl\n"
                "Or export your workbook as .csv and pass the .csv file instead."
            ) from exc

    @staticmethod
    def _read_excel_sheet_with_detected_header(path: Path, sheet_name: str) -> pd.DataFrame:
        preview = pd.read_excel(path, sheet_name=sheet_name, header=None, nrows=15, engine="openpyxl")
        header_row_index = Inventory._find_header_row_index(preview)
        if header_row_index is None:
            raise ValueError(f"Could not detect header row in sheet '{sheet_name}'")
        return pd.read_excel(path, sheet_name=sheet_name, header=header_row_index, engine="openpyxl")

    @staticmethod
    def _find_header_row_index(preview_df: pd.DataFrame) -> Optional[int]:
        required_markers = {
            "location",
            "code",
            "level",
            "contains (part/construct)",
            "part type",
            "plasmid conc. (ng/ul)",
            "length (bp)",
        }

        for idx in range(len(preview_df.index)):
            values = {
                str(value).strip().lower()
                for value in preview_df.iloc[idx].tolist()
                if pd.notna(value)
            }
            if required_markers.issubset(values):
                return idx
        return None

    @staticmethod
    def available_inventory_sheet_names(path: str | Path) -> List[str]:
        path = Path(path)
        if path.suffix.lower() == ".csv":
            return [path.name]

        try:
            excel_file = pd.ExcelFile(path, engine="openpyxl")
        except Exception as exc:
            # Inspect file header to decide how to proceed.
            try:
                with open(path, "rb") as fh:
                    magic = fh.read(4)
            except Exception:
                raise exc

            # ZIP-based .xlsx files start with 'PK\x03\x04'
            if magic.startswith(b"PK"):
                try:
                    excel_file = pd.ExcelFile(path, engine="openpyxl")
                except ImportError as iexc:
                    raise RuntimeError(
                        "Reading .xlsx files requires the 'openpyxl' package. "
                        "Install it with: pip install openpyxl\n"
                        "Or upload/export your inventory as .csv instead."
                    ) from iexc
                except Exception:
                    raise RuntimeError(
                        "The file appears to be an .xlsx workbook but could not be opened. "
                        "It may be corrupted or use an unsupported format. Try resaving as .xlsx or .csv."
                    ) from exc

            # Older Excel .xls files are compound file binary with magic 0xD0CF11E0
            elif magic == b"\xD0\xCF\x11\xE0":
                raise RuntimeError(
                    "The file appears to be an old Excel .xls file. "
                    "Pandas may require 'xlrd' to read this format, or resave as .xlsx/.csv."
                )

            else:
                # Not a binary Excel format — try parsing as CSV/text
                try:
                    _ = pd.read_csv(path)
                    return [path.name]
                except Exception:
                    raise exc
        valid_sheet_names: List[str] = []
        for name in excel_file.sheet_names:
            if name.lower() in {"general info", "sheet4"}:
                continue
            try:
                preview = pd.read_excel(path, sheet_name=name, header=None, nrows=15, engine="openpyxl")
                header_row_index = Inventory._find_header_row_index(preview)
                if header_row_index is not None:
                    valid_sheet_names.append(name)
            except Exception:
                continue
        return valid_sheet_names

    @staticmethod
    def _normalise_dataframe(df: pd.DataFrame) -> pd.DataFrame:
        present_columns_lower = {str(col).strip().lower() for col in df.columns}

        missing_columns = []
        for column in file_columns:
            if column.lower() not in present_columns_lower:
                missing_columns.append(column)
        if missing_columns:
            missing_text = "\n - ".join([""] + missing_columns)
            raise ValueError(f"Missing required columns:{missing_text}")

        # Flexible renaming (case-insensitive)
        rename_map = {}
        for col in df.columns:
            col_lower = str(col).strip().lower()
            for key, val in file_columns.items():
                if col_lower == key.lower():
                    rename_map[col] = val
        df = df.rename(columns=rename_map).copy()

        string_columns = [
            "code",
            "location",
            "level",
            "backbone",
            "ori",
            "resistance",
            "contains",
            "part_type",
            "a260_230",
            "a260_280",
            "sequenced",
            "five_prime_overhang",
            "three_prime_overhang",
            "description",
            "depositor",
            "benchling_file",
            "notes",
        ]
        for column in string_columns:
            if column in df.columns:
                df[column] = df[column].fillna("").astype(str).str.strip()

        numeric_columns = ["plasmid_concentration", "length_bp"]
        for column in numeric_columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

        # Ignore empty spreadsheet rows. Missing measurements are retained as NaN
        # and rejected later if a user selects that stock for a reaction.
        df = df[df["code"] != ""].copy()

        if "inventory_sheet" not in df.columns:
            df["inventory_sheet"] = ""

        level_map = {"0": "L0", "0.0": "L0", "1": "L1", "1.0": "L1", "2": "L2", "2.0": "L2"}
        df["level"] = df["level"].map(
            lambda value: level_map.get(str(value).strip().casefold(), str(value).strip().upper())
        )

        return df

    @staticmethod
    def _rows_to_items(df: pd.DataFrame) -> List[DNAItem]:
        items: List[DNAItem] = []

        for _, row in df.iterrows():
            items.append(
                DNAItem(
                    code=row["code"],
                    location=row.get("location", ""),
                    level=row.get("level", ""),
                    backbone=row.get("backbone", ""),
                    resistance=row.get("resistance", ""),
                    ori=row.get("ori", ""),
                    contains=row.get("contains", ""),
                    part_type=row.get("part_type", ""),
                    length_bp=float(row["length_bp"]),
                    concentration_ng_ul=float(row["plasmid_concentration"]),
                    a260_230=row.get("a260_230", ""),
                    a260_280=row.get("a260_280", ""),
                    sequenced=row.get("sequenced", ""),
                    five_prime_overhang=row.get("five_prime_overhang", ""),
                    three_prime_overhang=row.get("three_prime_overhang", ""),
                    notes=row.get("notes", ""),
                    inventory_sheet=row.get("inventory_sheet", ""),
                )
            )

        return items

    def add_item(self, item: DNAItem) -> None:
        if any(
            existing.code.strip().casefold() == item.code.strip().casefold()
            and existing.location.strip().casefold() == item.location.strip().casefold()
            and existing.inventory_sheet.strip().casefold() == item.inventory_sheet.strip().casefold()
            for existing in self.items
        ):
            raise ValueError(
                f"Duplicate inventory code/location record: {item.code} "
                f"[{item.inventory_sheet} / {item.location}]."
            )
        self.items.append(item)

    def find_by_code(self, code: str) -> Optional[DNAItem]:
        matches = self.find_all_by_code(code)
        return matches[0] if len(matches) == 1 else None

    def find_all_by_code(self, code: str) -> List[DNAItem]:
        code = code.strip().casefold()
        return [item for item in self.items if item.code.strip().casefold() == code]

    def search(self, text: str) -> List[DNAItem]:
        text = text.strip().lower()
        matches: List[DNAItem] = []
        for item in self.items:
            haystack = " ".join(
                [
                    item.code,
                    item.contains,
                    item.part_type,
                    item.level,
                    item.backbone,
                    item.resistance,
                    item.location,
                ]
            ).lower()
            if text in haystack:
                matches.append(item)
        return matches

    def list_codes(self, limit: int = 20) -> None:
        for item in self.items[:limit]:
            print(f" - {item.label}")


# Inserted choose_item_by_code function

def create_manual_dna_item(
    code: str,
    contains: str,
    part_type: str,
    length_bp: float,
    concentration_ng_ul: float,
    level: str = "custom",
    location: str = "manual",
    backbone: str = "",
    resistance: str = "",
    notes: str = "",
    five_prime_overhang: str = "",
    three_prime_overhang: str = "",
) -> DNAItem:
    return CustomDNAInput(
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
    ).to_dna_item()

def ng_for_fmol(length_bp: float, target_fmol: float) -> float:
    if not math.isfinite(length_bp) or length_bp <= 0:
        raise ValueError("length_bp must be a positive finite value in bp")
    if not math.isfinite(target_fmol) or target_fmol <= 0:
        raise ValueError("target_fmol must be a positive finite value")
    return target_fmol * length_bp * 0.00065


def volume_for_fmol(length_bp: float, concentration_ng_ul: float, target_fmol: float) -> float:
    if not math.isfinite(concentration_ng_ul) or concentration_ng_ul <= 0:
        raise ValueError("concentration_ng_ul must be a positive finite value in ng/uL")
    return ng_for_fmol(length_bp, target_fmol) / concentration_ng_ul


def build_reaction(
    backbone: DNAItem,
    inserts: List[DNAItem],
    backbone_fmol: float,
    insert_ratio: float = 1.0,
) -> List[ReactionComponent]:
    if not math.isfinite(backbone_fmol) or backbone_fmol <= 0:
        raise ValueError("backbone_fmol must be a positive finite value")
    if not math.isfinite(insert_ratio) or insert_ratio <= 0:
        raise ValueError("insert_ratio must be a positive finite value")
    components = [ReactionComponent(item=backbone, target_fmol=backbone_fmol)]
    insert_fmol = backbone_fmol * insert_ratio

    for insert in inserts:
        components.append(ReactionComponent(item=insert, target_fmol=insert_fmol))

    return components


def reaction_to_dataframe(components: List[ReactionComponent]) -> pd.DataFrame:
    rows = []
    for component in components:
        rows.append(
            {
                "Code": component.item.code,
                "Contains": component.item.contains,
                "Part type": component.item.part_type,
                "Role": component.role,
                "Level": component.item.level,
                "Length (bp)": round(component.item.length_bp, 1),
                "Conc. (ng/uL)": round(component.item.concentration_ng_ul, 2),
                "Inventory sheet": component.item.inventory_sheet,
                "Target (fmol)": round(component.target_fmol, 2),
                "Required DNA (ng)": round(component.required_ng, 3),
                "Volume to add (uL)": round(component.required_volume_ul, 3),
                "Location": component.item.location,
            }
        )
    return pd.DataFrame(rows)


def calculate_mq_water(components: List[ReactionComponent], settings: ReactionProtocol) -> float:
    total_dna_volume = sum(component.required_volume_ul for component in components)
    water_ul = settings.total_reaction_ul - settings.fixed_volume_ul - total_dna_volume
    if water_ul < -1e-9:
        raise ValueError(
            f"Reaction volumes do not fit: DNA and fixed reagents require "
            f"{total_dna_volume + settings.fixed_volume_ul:.3f} uL, exceeding "
            f"the {settings.total_reaction_ul:.3f} uL total volume."
        )
    return max(0.0, water_ul)



def build_level0_reaction_rows(
    vector: DNAItem,
    inserts: List[DNAItem],
    fmol_target: float,
    settings: ReactionProtocol,
) -> pd.DataFrame:
    rows = []

    for index, insert in enumerate(inserts, start=1):
        insert_component = ReactionComponent(item=insert, target_fmol=fmol_target)
        vector_component = ReactionComponent(item=vector, target_fmol=fmol_target)
        mq_water_ul = calculate_mq_water(
            [vector_component, insert_component], settings
        )

        rows.append(
            {
                "Insert Tube label": "",
                "Sample Number": index,
                "Sample Name": f"{insert.contains}_{vector.code}",
                "Insert Name": insert.contains,
                "Insert Length (bp)": round(insert.length_bp, 1),
                "Insert Conc (ng/µL)": round(insert.concentration_ng_ul, 2),
                "Vector Name": vector.code,
                "Vector Length (bp)": round(vector.length_bp, 1),
                "Vector Conc (ng/µL)": round(vector.concentration_ng_ul, 2),
                "fmol Target": round(fmol_target, 2),
                "Insert Vol (µL)": round(insert_component.required_volume_ul, 3),
                "Vector Vol (µL)": round(vector_component.required_volume_ul, 3),
                "MQ water (µL)": round(mq_water_ul, 3),
                "Completed": "",
            }
        )

    return pd.DataFrame(rows)

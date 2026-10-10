from pathlib import Path
import math

import pandas as pd
import pytest

from moclo_core import (
    DNAItem,
    Inventory,
    ReactionComponent,
    ReactionProtocol,
    calculate_mq_water,
    file_columns,
    volume_for_fmol,
)


def dna_item(code="p1", length_bp=1000, concentration_ng_ul=2.0):
    return DNAItem(
        code=code,
        location="A1",
        level="L1",
        backbone="",
        resistance="",
        contains=code,
        part_type="CDS",
        length_bp=length_bp,
        concentration_ng_ul=concentration_ng_ul,
    )


def test_fmol_to_mass_and_volume_uses_ds_dna_approximation():
    component = ReactionComponent(dna_item(), target_fmol=20)

    assert component.required_ng == pytest.approx(13.0)
    assert component.required_volume_ul == pytest.approx(6.5)
    assert volume_for_fmol(1000, 2.0, 20) == pytest.approx(6.5)


def test_missing_or_nonpositive_concentration_is_rejected():
    component = ReactionComponent(dna_item(concentration_ng_ul=0), target_fmol=20)

    with pytest.raises(ValueError, match="positive finite value in ng/uL"):
        _ = component.required_volume_ul


def test_reaction_with_impossible_total_volume_is_rejected():
    vector = ReactionComponent(dna_item("vector", concentration_ng_ul=1), target_fmol=20)
    insert = ReactionComponent(dna_item("insert", concentration_ng_ul=1), target_fmol=40)
    protocol = ReactionProtocol(
        total_reaction_ul=10,
        buffer_ul=1,
        ligase_ul=0.5,
        enzyme_ul=0.5,
    )

    with pytest.raises(ValueError, match="exceeding the 10.000 uL total volume"):
        calculate_mq_water([vector, insert], protocol)


def test_duplicate_inventory_code_and_location_are_rejected_case_insensitively():
    with pytest.raises(ValueError, match="Duplicate inventory code"):
        Inventory([dna_item("pA"), dna_item("pa")])


def test_same_code_at_different_stock_locations_is_allowed():
    first = dna_item("pA")
    second = dna_item("pA")
    second.location = "B2"

    inventory = Inventory([first, second])

    assert len(inventory.items) == 2
    assert inventory.find_by_code("pA") is None
    assert len(inventory.find_all_by_code("pA")) == 2


def test_same_code_and_location_on_different_sheets_is_allowed():
    first = dna_item("pA")
    second = dna_item("pA")
    second.inventory_sheet = "Second plate"

    inventory = Inventory([first, second])

    assert len(inventory.items) == 2


def complete_row():
    row = {header: "" for header in file_columns}
    row.update(
        {
            "Location": "A1",
            "Code": "vector-1",
            "Level": "L1",
            "Plasmid backbone": "DVK",
            "Contains (part/construct)": "synthetic vector",
            "Part type": "Destination vector",
            "5' flank": "A",
            "3' flank": "E",
            "Plasmid conc. (ng/uL)": 10,
            "Length (bp)": 5000,
        }
    )
    return row


def test_excel_loader_detects_offset_header_and_case_insensitive_columns(tmp_path: Path):
    columns = list(file_columns)
    reversed_columns = columns[::-1]
    frame = pd.DataFrame([complete_row()], columns=reversed_columns)
    frame.columns = [column.swapcase() for column in reversed_columns]
    path = tmp_path / "formatted_inventory.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        frame.to_excel(writer, sheet_name="MoClo Dest Vectors", index=False, startrow=2)
        writer.sheets["MoClo Dest Vectors"]["A1"] = "Synthetic inventory title"

    sheet_names = Inventory.available_inventory_sheet_names(path)
    loaded = Inventory.from_excel_sheets(path, sheet_names)

    assert sheet_names == ["MoClo Dest Vectors"]
    assert len(loaded.items) == 1
    assert loaded.items[0].five_prime_overhang == "A"
    assert loaded.items[0].three_prime_overhang == "E"


def test_excel_loader_retains_missing_concentration_for_use_time_validation(tmp_path: Path):
    row = complete_row()
    row["Plasmid conc. (ng/uL)"] = None
    path = tmp_path / "missing_concentration.xlsx"
    pd.DataFrame([row], columns=list(file_columns)).to_excel(
        path, sheet_name="MoClo Dest Vectors", index=False
    )

    loaded = Inventory.from_excel_sheets(path, ["MoClo Dest Vectors"])

    assert len(loaded.items) == 1
    assert math.isnan(loaded.items[0].concentration_ng_ul)
    with pytest.raises(ValueError, match="positive finite value in ng/uL"):
        _ = ReactionComponent(loaded.items[0], target_fmol=20).required_volume_ul


def test_excel_loader_treats_malformed_measurements_as_missing(tmp_path: Path):
    row = complete_row()
    row["Plasmid conc. (ng/uL)"] = "115..2"
    row["Length (bp)"] = None
    path = tmp_path / "malformed_measurements.xlsx"
    pd.DataFrame([row], columns=list(file_columns)).to_excel(
        path, sheet_name="MoClo Dest Vectors", index=False
    )

    loaded = Inventory.from_excel_sheets(path, ["MoClo Dest Vectors"])

    assert math.isnan(loaded.items[0].concentration_ng_ul)
    assert math.isnan(loaded.items[0].length_bp)
    with pytest.raises(ValueError, match="Invalid DNA length for vector-1"):
        _ = ReactionComponent(loaded.items[0], target_fmol=20).required_ng


def test_excel_loader_normalizes_numeric_assembly_levels(tmp_path: Path):
    row = complete_row()
    row["Level"] = 0
    path = tmp_path / "numeric_level.xlsx"
    pd.DataFrame([row], columns=list(file_columns)).to_excel(
        path, sheet_name="MoClo Dest Vectors", index=False
    )

    loaded = Inventory.from_excel_sheets(path, ["MoClo Dest Vectors"])

    assert loaded.items[0].level == "L0"

"""Streamlit user interface for MoCloPlanner. Run with: streamlit run moclo_gui.py"""

import math
from itertools import product
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import streamlit as st

from moclo_core import (
    Inventory, ReactionProtocol, build_reaction, calculate_mq_water,
    reaction_to_dataframe,
)
from moclo_assembly import (
    clean, get_next_insert_index,
    get_valid_parts_for_position, part_identifier, dedupe_parts,
    chain_signature, construct_label,
)
from moclo_standards import (
    DEFAULT_STANDARD, STANDARDS, get_expected_flanks,
    validate_part_chain, validate_vector,
)

st.set_page_config(page_title="MoCloPlanner", layout="wide")

st.title("🧬 MoCloPlanner")
standard_name = st.selectbox(
    "Assembly standard",
    options=list(STANDARDS),
    index=list(STANDARDS).index(DEFAULT_STANDARD.name),
)
assembly_standard = STANDARDS[standard_name]

# -----------------------------
# File upload
# -----------------------------
uploaded_file = st.file_uploader("Upload inventory Excel file", type=["xlsx"])

if uploaded_file is None:
    st.stop()

# Save temp file
temp_suffix = Path(uploaded_file.name).suffix or ".xlsx"
temp_path = Path(f"temp_inventory{temp_suffix}")

# Write uploaded file safely
try:
    data = uploaded_file.getbuffer()
except Exception:
    uploaded_file.seek(0)
    data = uploaded_file.read()

with open(temp_path, "wb") as f:
    f.write(data)

# Sanity check
import os
if os.path.getsize(temp_path) == 0:
    st.error("Uploaded file is empty or failed to save correctly.")
    st.stop()

# -----------------------------
# Assembly level
# -----------------------------
assembly_level = st.selectbox(
    "Assembly level",
    ["L0", "L1", "L2"]
)

# initialise session state
if "selected_parts" not in st.session_state:
    st.session_state.selected_parts = []

if "manual_parts" not in st.session_state:
    st.session_state.manual_parts = []

required_part_level = assembly_standard.inventory_part_level(assembly_level)

# -----------------------------
# Manual part entry
# -----------------------------
with st.expander("Manually add part", expanded=False):
    st.caption("Adds a custom part into the current candidate pool for debugging or quick testing.")

    with st.form("manual_part_form", clear_on_submit=True):
        manual_label = st.text_input("Label")
        manual_part_type = st.text_input("Part type")
        manual_length = st.number_input("Length (bp)", min_value=0.0, value=0.0, step=1.0)
        manual_concentration = st.number_input("Plasmid conc. (ng/uL)", min_value=0.0, value=0.0, step=0.1)
        manual_five_prime = st.text_input("5' flank")
        manual_three_prime = st.text_input("3' flank")
        submitted_manual_part = st.form_submit_button("Add manual part")

    if submitted_manual_part:
        if not manual_label.strip():
            st.warning("Enter a label for the manual part.")
        elif not manual_part_type.strip():
            st.warning("Enter a part type for the manual part.")
        elif manual_length <= 0:
            st.warning("Enter a positive DNA length in bp.")
        elif manual_concentration <= 0:
            st.warning("Enter a positive DNA concentration in ng/uL.")
        else:
            manual_code = f"manual_{len(st.session_state.manual_parts) + 1}"
            st.session_state.manual_parts.append(
                SimpleNamespace(
                    label=manual_label.strip(),
                    code=manual_code,
                    level=required_part_level,
                    part_type=manual_part_type.strip(),
                    length_bp=manual_length,
                    concentration_ng_ul=manual_concentration,
                    ori="",
                    contains=manual_label.strip(),
                    five_prime_overhang=manual_five_prime.strip().upper(),
                    three_prime_overhang=manual_three_prime.strip().upper(),
                    a260_230="",
                    a260_280="",
                    sequenced="",
                )
            )
            st.success(f"Added manual part: {manual_label.strip()}")

if st.session_state.manual_parts:
    st.write("Manual parts loaded:")
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "label": p.label,
                    "part_type": p.part_type,
                    "concentration_ng_ul": p.concentration_ng_ul,
                    "level": p.level,
                }
                for p in st.session_state.manual_parts
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )
# -----------------------------
# Load inventory
# -----------------------------
available_sheets = Inventory.available_inventory_sheet_names(temp_path)

selected_sheets = st.multiselect(
    "Select inventory sheets",
    options=available_sheets,
    default=available_sheets
)

if not selected_sheets:
    st.warning("Select at least one sheet")
    st.stop()

try:
    inventory = Inventory.from_excel_sheets(temp_path, selected_sheets)
except ValueError as exc:
    st.error(f"Inventory validation failed: {exc}")
    st.stop()

st.success(f"Loaded {len(inventory.items)} items")

incomplete_stock_rows = []
for item in inventory.items:
    missing_fields = []
    for field_name, value, unit in (
        ("length", item.length_bp, "bp"),
        ("concentration", item.concentration_ng_ul, "ng/uL"),
    ):
        try:
            is_valid = math.isfinite(float(value)) and float(value) > 0
        except (TypeError, ValueError):
            is_valid = False
        if not is_valid:
            missing_fields.append(f"{field_name} ({unit})")
    if missing_fields:
        incomplete_stock_rows.append(
            {
                "Code": item.code,
                "Inventory sheet": item.inventory_sheet,
                "Location": item.location,
                "Missing or invalid": ", ".join(missing_fields),
            }
        )

if incomplete_stock_rows:
    st.warning(
        f"{len(incomplete_stock_rows)} stock row(s) lack a positive length or concentration. "
        "They can be reviewed below, but reaction calculation will reject a selected row "
        "until its measurements are corrected."
    )
    st.dataframe(pd.DataFrame(incomplete_stock_rows), use_container_width=True, hide_index=True)

# -----------------------------
# Backbone selection
# -----------------------------
# Load "MoClo Dest Vectors" sheet specifically for backbone selection
backbone_sheets = ["MoClo Dest Vectors"] if "MoClo Dest Vectors" in available_sheets else selected_sheets
try:
    backbone_inventory = (
        Inventory.from_excel_sheets(temp_path, backbone_sheets)
        if backbone_sheets else inventory
    )
except ValueError as exc:
    st.error(f"Backbone inventory validation failed: {exc}")
    st.stop()

all_items = inventory.items + st.session_state.manual_parts
backbone_candidates = backbone_inventory.items


st.write("Assembly level:", assembly_level)
st.write("Allowed backbone families:", assembly_standard.allowed_vector_families(assembly_level))
st.write(
    pd.DataFrame([
        {
            "label": v.label,
            "part_type": getattr(v, "part_type", ""),
            "plasmid_backbone": getattr(v, "backbone", ""),
        }
        for v in backbone_candidates
    ])
)
vectors = []
invalid_vector_messages = []
for candidate in backbone_candidates:
    part_type = clean(getattr(candidate, "part_type", ""))
    if "VECTOR" not in part_type and "BACKBONE" not in part_type:
        continue
    if clean(getattr(candidate, "backbone", "")) not in assembly_standard.allowed_vector_families(assembly_level):
        continue
    try:
        validate_vector(candidate, assembly_level, assembly_standard)
    except ValueError as exc:
        invalid_vector_messages.append(str(exc))
        continue
    vectors.append(candidate)
if not vectors:
    st.warning("No valid vectors found for this assembly level.")
    if invalid_vector_messages:
        st.caption("Check explicit 5' and 3' flank values in the vector rows.")
        for message in invalid_vector_messages:
            st.write(f"- {message}")
    st.stop()
if invalid_vector_messages:
    st.warning(f"Ignored {len(invalid_vector_messages)} vector row(s) with invalid {assembly_standard.name} settings.")

backbone_selection = st.selectbox(
    "Select backbone",
    options=vectors,
    format_func=lambda x: x.label
)
if backbone_selection is None:
    st.error("Select a valid backbone before planning an assembly.")
    st.stop()
# `selectbox` is typed as possibly returning None, while the earlier empty-list
# guard guarantees that at least one backbone option is available.
assert backbone_selection is not None
backbone = backbone_selection

# -----------------------------
# Insert selection (CIDAR overhang-aware)
# -----------------------------

# reset if the selected backbone stock changes
backbone_selection_key = part_identifier(backbone)
if "last_backbone" not in st.session_state or st.session_state.last_backbone != backbone_selection_key:
    st.session_state.selected_parts = []
    st.session_state.last_backbone = backbone_selection_key


parts = [
     i for i in all_items
     if hasattr(i, "level")
     and str(i.level).strip().upper() == required_part_level
 ]

if not parts:
    parts = inventory.items

st.subheader(f"Stepwise part selection ({assembly_standard.name} order)")

# ============================================================
# LEVEL 0: Manual insert entry (PCR products / donor plasmids)
# ============================================================
if assembly_standard.uses_manual_inserts(assembly_level):
    st.info("Level 0 assemblies require manually entered insert parts (linear PCR products or donor plasmids).")
    
    # show current chain
    if st.session_state.selected_parts:
        st.write("Current inserts:")
        for idx, p in enumerate(st.session_state.selected_parts, start=1):
            if p is None:
                st.write(f"{idx}. [missing]")
            else:
                st.write(f"{idx}. {p.label} ({p.part_type})")

        remove_col1, remove_col2 = st.columns(2)

        removable_indexes = [
            i for i, p in enumerate(st.session_state.selected_parts)
            if p is not None
        ]

        with remove_col1:
            if removable_indexes:
                part_to_remove = st.selectbox(
                    "Remove an insert",
                    options=removable_indexes,
                    format_func=lambda i: st.session_state.selected_parts[i].label,
                    key="remove_part_select_l0",
                )
                if st.button("Remove selected insert"):
                    st.session_state.selected_parts[part_to_remove] = None
                    st.rerun()

        with remove_col2:
            if st.button("Remove last insert"):
                while st.session_state.selected_parts and st.session_state.selected_parts[-1] is None:
                    st.session_state.selected_parts.pop()
                if st.session_state.selected_parts:
                    st.session_state.selected_parts.pop()
                st.rerun()

    # Manual L0 insert entry
    with st.expander("Add L0 insert (PCR product or donor plasmid)", expanded=True):
        with st.form("l0_insert_form", clear_on_submit=True):
            insert_label = st.text_input("Insert label/name")
            insert_type = st.selectbox("Insert type", ["PCR product", "Donor plasmid", "Other"])
            insert_length = st.number_input("Length (bp)", min_value=0.0, value=0.0, step=1.0)
            insert_concentration = st.number_input("Concentration (ng/uL)", min_value=0.0, value=0.0, step=0.1)
            insert_five_prime = st.text_input("5' overhang", placeholder="e.g., A")
            insert_three_prime = st.text_input("3' overhang", placeholder="e.g., B")
            submitted_l0_insert = st.form_submit_button("Add insert")

        if submitted_l0_insert:
            if not insert_label.strip():
                st.warning("Enter a label for the insert.")
            elif not insert_five_prime.strip() or not insert_three_prime.strip():
                st.warning("Enter both 5' and 3' overhangs.")
            elif insert_length <= 0:
                st.warning("Enter a positive DNA length in bp.")
            elif insert_concentration <= 0:
                st.warning("Enter a positive DNA concentration in ng/uL.")
            else:
                insert_code = f"L0_insert_{len(st.session_state.selected_parts) + 1}"
                st.session_state.selected_parts.append(
                    SimpleNamespace(
                        label=insert_label.strip(),
                        code=insert_code,
                        level="L0",
                        part_type=insert_type,
                        length_bp=insert_length,
                        concentration_ng_ul=insert_concentration,
                        ori="",
                        contains=insert_label.strip(),
                        five_prime_overhang=insert_five_prime.strip().upper(),
                        three_prime_overhang=insert_three_prime.strip().upper(),
                        a260_230="",
                        a260_280="",
                        sequenced="",
                    )
                )
                st.success(f"Added L0 insert: {insert_label.strip()}")

# ============================================================
# LEVEL 1 & 2: Inventory-based part selection (CIDAR order)
# ============================================================
else:
    st.write("Current chain:")
    for idx, p in enumerate(st.session_state.selected_parts, start=1):
        if p is None:
            st.write(f"{idx}. [missing]")
        else:
            st.write(f"{idx}. {p.label}")

    remove_col1, remove_col2 = st.columns(2)

    removable_indexes = [
        i for i, p in enumerate(st.session_state.selected_parts)
        if p is not None
    ]

    with remove_col1:
        if removable_indexes:
            part_to_remove = st.selectbox(
                "Remove a part",
                options=removable_indexes,
                format_func=lambda i: st.session_state.selected_parts[i].label,
                key="remove_part_select",
            )
            if st.button("Remove selected part"):
                # Keep downstream entries and mark this slot missing for replacement.
                st.session_state.selected_parts[part_to_remove] = None
                st.rerun()

    with remove_col2:
        if st.button("Remove last part"):
            while st.session_state.selected_parts and st.session_state.selected_parts[-1] is None:
                st.session_state.selected_parts.pop()
            if st.session_state.selected_parts:
                st.session_state.selected_parts.pop()
            st.rerun()

    # Get expected flanks for the next insert
    next_insert_index = get_next_insert_index(st.session_state.selected_parts)
    expected_5, expected_3 = get_expected_flanks(
        backbone,
        next_insert_index,
        assembly_level,
        assembly_standard,
    )

    if expected_5 is None or expected_3 is None:
        st.warning("Could not determine valid flanks from backbone")

    else:
        st.write(f"Next required flanks: **5' = {expected_5}, 3' = {expected_3}**")
        show_all_parts = st.checkbox(
            "Show all parts (bypass expected filtering)",
            value=False,
            key="show_all_parts",
        )
        _, _, valid_parts = get_valid_parts_for_position(
            parts,
            backbone,
            next_insert_index,
            assembly_level,
            show_all_parts=show_all_parts,
            standard=assembly_standard,
        )

        if show_all_parts:
            st.info("Showing every loaded part; final flank validation still applies before calculation.")

        if not valid_parts:
            st.warning(f"No parts found with 5' = {expected_5} and 3' = {expected_3}")

        else:
            def na(value):
                if pd.isna(value):
                    return "Missing"
                text = str(value).strip()
                return text if text else "Missing"

            valid_parts_df = pd.DataFrame(
                {
                    "Part": [p.label for p in valid_parts],
                    "Ori": [na(getattr(p, "ori", "")) for p in valid_parts],
                    "Plasmid conc. (ng/uL)": [na(getattr(p, "concentration_ng_ul", "")) for p in valid_parts],
                    "A260/230": [na(getattr(p, "a260_230", "")) for p in valid_parts],
                    "A260/280": [na(getattr(p, "a260_280", "")) for p in valid_parts],
                    "Sequenced": [na(getattr(p, "sequenced", "")) for p in valid_parts],
                }
            )
            st.dataframe(valid_parts_df, use_container_width=True, hide_index=True)

            valid_part_labels = [p.label for p in valid_parts]

            with st.form(key=f"add_part_form_{next_insert_index}"):
                selected_part_index = st.selectbox(
                    "Select next part",
                    options=list(range(len(valid_parts))),
                    format_func=lambda idx: valid_part_labels[idx],
                    key=f"part_select_{next_insert_index}",
                )
                add_part_clicked = st.form_submit_button("Add part")

            selected = valid_parts[selected_part_index]

            if add_part_clicked:
                if next_insert_index < len(st.session_state.selected_parts):
                    st.session_state.selected_parts[next_insert_index] = selected
                else:
                    st.session_state.selected_parts.append(selected)

# final selected parts list
selected_parts = st.session_state.selected_parts
selected_part_groups = [[part] for part in selected_parts if part is not None]
batch_variant_mode = False
planned_assemblies = 1 if selected_part_groups else 0

if (
    assembly_standard.supports_variant_assemblies(assembly_level)
    and selected_parts
    and not any(part is None for part in selected_parts)
):
    st.subheader("Variant assemblies")
    batch_variant_mode = st.checkbox(
        "Optional: calculate multiple assemblies from chain-position variants",
        value=False,
        key="batch_variant_mode",
    )

    if batch_variant_mode:
        st.caption(
            "Select one or more allowed parts for each chain position. "
            "The app will calculate every resulting combination with the same backbone."
        )

        selected_part_groups = []
        planned_assemblies = 1
        variant_key_prefix = chain_signature(backbone, selected_parts)

        for position_index, current_part in enumerate(selected_parts):
            expected_5, expected_3, candidate_parts = get_valid_parts_for_position(
                parts,
                backbone,
                position_index,
                assembly_level,
                standard=assembly_standard,
            )
            position_options = dedupe_parts([current_part] + candidate_parts)
            option_labels = [part.label for part in position_options]
            default_indices = [
                index for index, option in enumerate(position_options)
                if part_identifier(option) == part_identifier(current_part)
            ]
            position_title = f"Position {position_index + 1} ({expected_5}->{expected_3}) variants"

            selected_indices = st.multiselect(
                position_title,
                options=list(range(len(position_options))),
                default=default_indices,
                format_func=lambda index, labels=option_labels: labels[index],
                key=f"variant_position_{position_index}_{variant_key_prefix}",
            )

            selected_options = [position_options[index] for index in selected_indices]
            selected_part_groups.append(selected_options)

            if not selected_options:
                planned_assemblies = 0
                st.warning(f"Select at least one part for position {position_index + 1}.")
            else:
                planned_assemblies *= len(selected_options)

        if planned_assemblies:
            st.write(f"Planned assemblies: **{planned_assemblies}**")
            if planned_assemblies > 24:
                st.warning(
                    "This will generate many reaction plans. "
                    "Keep the variant count small unless you really want the full combination set."
                )

# -----------------------------
# Reaction settings
# -----------------------------
col1, col2, col3 = st.columns(3)

with col1:
    backbone_fmol = st.number_input("Backbone target (fmol)", min_value=0.01, value=20.0)

with col2:
    insert_ratio = st.number_input("Insert ratio (relative to backbone)", min_value=0.01, value=2.0)

with col3:
    total_volume = st.number_input("Total reaction volume (uL)", min_value=0.1, value=10.0)

with st.expander("Reaction composition", expanded=False):
    st.caption("Set these for your laboratory protocol. They are reaction settings, separate from CIDAR assembly rules.")
    buffer_volume = st.number_input("Buffer volume (uL)", min_value=0.0, value=1.0, step=0.1)
    ligase_volume = st.number_input("Ligase volume (uL)", min_value=0.0, value=0.5, step=0.1)
    enzyme_name = st.text_input("Restriction enzyme name", value="Restriction enzyme")
    enzyme_volume = st.number_input("Restriction enzyme volume (uL)", min_value=0.0, value=0.5, step=0.1)

try:
    settings = ReactionProtocol(
        total_reaction_ul=total_volume,
        buffer_ul=buffer_volume,
        ligase_ul=ligase_volume,
        enzyme_ul=enzyme_volume,
        backbone_fmol=backbone_fmol,
        insert_ratio=insert_ratio,
        enzyme_name=enzyme_name,
    )
except ValueError as exc:
    st.error(f"Invalid reaction settings: {exc}")
    st.stop()

# -----------------------------
# Run calculation
# -----------------------------
if st.button("Calculate reaction"):
    if not selected_parts:
        st.warning("Select at least one insert")
        st.stop()

    if any(p is None for p in selected_parts):
        st.warning("Fill all missing parts in the chain before calculating")
        st.stop()

    try:
        validate_part_chain(backbone, selected_parts, assembly_level, assembly_standard)
    except ValueError as exc:
        st.error(f"Assembly validation failed: {exc}")
        st.stop()

    if batch_variant_mode:
        if not selected_part_groups or any(not group for group in selected_part_groups):
            st.warning("Select at least one part in every variant position before calculating.")
            st.stop()

        variant_combinations = list(product(*selected_part_groups))
        summary_rows = []
        detailed_frames = []

        for assembly_index, inserts in enumerate(variant_combinations, start=1):
            try:
                validate_part_chain(backbone, inserts, assembly_level, assembly_standard)
                components = build_reaction(
                    backbone=backbone,
                    inserts=list(inserts),
                    backbone_fmol=settings.backbone_fmol,
                    insert_ratio=settings.insert_ratio,
                )
                reaction_df = reaction_to_dataframe(components)
                mq_water = calculate_mq_water(components, settings)
            except ValueError as exc:
                st.error(f"Assembly {assembly_index} could not be planned: {exc}")
                st.stop()
            reaction_df.insert(0, "Assembly", f"Assembly {assembly_index}")
            reaction_df.insert(1, "Construct", construct_label(backbone, inserts))
            detailed_frames.append(reaction_df)

            total_dna = reaction_df["Volume to add (uL)"].sum()

            summary_row = {
                "Assembly": f"Assembly {assembly_index}",
                "Backbone": backbone.code,
                "Construct": construct_label(backbone, inserts),
                "Total DNA volume (uL)": round(total_dna, 3),
                "MQ water (uL)": round(mq_water, 3),
            }
            for position_index, insert in enumerate(inserts, start=1):
                expected_5, expected_3 = get_expected_flanks(
                    backbone, position_index - 1, assembly_level, assembly_standard
                )
                summary_row[f"Insert {position_index} ({expected_5}->{expected_3})"] = insert.code
            summary_rows.append(summary_row)

        st.subheader("Assembly summary")
        st.dataframe(pd.DataFrame(summary_rows), use_container_width=True, hide_index=True)

        st.subheader("Reaction plans")
        st.dataframe(pd.concat(detailed_frames, ignore_index=True), use_container_width=True)
        st.write(f"**Assemblies calculated:** {len(variant_combinations)}")
        st.caption(
            f"Protocol composition: {settings.total_reaction_ul:g} uL total; "
            f"{settings.buffer_ul:g} uL buffer; {settings.ligase_ul:g} uL ligase; "
            f"{settings.enzyme_ul:g} uL {settings.enzyme_name}; target "
            f"{settings.backbone_fmol:g} fmol backbone, insert ratio {settings.insert_ratio:g}."
        )
        st.stop()

    try:
        components = build_reaction(
            backbone=backbone,
            inserts=selected_parts,
            backbone_fmol=settings.backbone_fmol,
            insert_ratio=settings.insert_ratio,
        )
        df = reaction_to_dataframe(components)
        water_volume = calculate_mq_water(components, settings)
    except ValueError as exc:
        st.error(f"Reaction could not be planned: {exc}")
        st.stop()

    st.subheader("Reaction plan")
    st.dataframe(df, use_container_width=True)

    total_dna = df["Volume to add (uL)"].sum()
    st.write(f"**Total DNA volume:** {total_dna:.3f} µL")
    st.write(f"**Water to add:** {water_volume:.3f} µL")
    st.caption(
        f"Protocol composition: {settings.total_reaction_ul:g} uL total; "
        f"{settings.buffer_ul:g} uL buffer; {settings.ligase_ul:g} uL ligase; "
        f"{settings.enzyme_ul:g} uL {settings.enzyme_name}; target "
        f"{settings.backbone_fmol:g} fmol backbone, insert ratio {settings.insert_ratio:g}."
    )

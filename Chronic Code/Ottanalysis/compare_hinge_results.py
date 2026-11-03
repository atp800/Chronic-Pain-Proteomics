"""Compare completed hinge runs using the current application configuration.

Call compare_hinge_results(result_files, CONFIG) from main.py after the hinge
loop. result_files maps each successfully completed group to its result path.
This prevents a failed current run from being replaced by an old workbook.

With two groups, Common_Significant is their significant overlap. With more
groups, it is the intersection across every selected group. A pooled run or
one selected group has no between-group overlap to calculate.
"""

from pathlib import Path

import pandas as pd

from config import CONFIG


def _read_hinge_results(path):
    """Read the standard result schema written by hinge_model.py."""
    results = pd.read_excel(
        path,
        sheet_name="All_Model_Results",
        dtype={"Protein": "string"},
    )
    required = {"Protein", "Slope_2_Change_PValue"}
    missing = required - set(results.columns)
    if missing:
        raise ValueError(f"{path}: missing result columns {sorted(missing)}")

    results["Protein"] = results["Protein"].str.strip()
    if results["Protein"].fillna("").eq("").any():
        raise ValueError(f"{path}: blank protein identifiers found")
    if results["Protein"].duplicated().any():
        raise ValueError(f"{path}: duplicate protein identifiers found")

    p_values = pd.to_numeric(
        results["Slope_2_Change_PValue"], errors="raise"
    )
    if not p_values.dropna().between(0, 1).all():
        raise ValueError(f"{path}: p-values must be between 0 and 1")
    results["Slope_2_Change_PValue"] = p_values
    return results


def compare_hinge_results(result_files, config=None):
    """Save the significant overlap for all configured hinge groups.

    Parameters
    ----------
    result_files : mapping
        Group name -> workbook path, collected from successful calls to
        run_hinge_model during this run. Do not populate it by scanning old
        output folders.
    config : mapping, optional
        The live CONFIG dictionary. Defaults to the imported CONFIG object.

    Returns
    -------
    pathlib.Path or None
        The saved summary path, or None when comparison was skipped.
    """
    if config is None:
        config = CONFIG

    groups = [
        str(group)
        for group in (config.get("HINGE_GROUPS_TO_RUN") or [])
        if group is not None
    ]
    if len(groups) != len(set(groups)):
        raise ValueError("HINGE_GROUPS_TO_RUN contains duplicate group names")
    if len(groups) < 2:
        print("Hinge comparison skipped: fewer than two separate groups selected.")
        return None

    threshold = float(config["P_THRESHOLD"])
    if not 0 <= threshold <= 1:
        raise ValueError("P_THRESHOLD must be between 0 and 1")
    if not config.get("OUTPUT_FOLDER"):
        raise ValueError("OUTPUT_FOLDER is required")

    paths = {str(group): Path(path) for group, path in result_files.items()}
    missing_groups = [
        group for group in groups
        if group not in paths or not paths[group].is_file()
    ]
    if missing_groups:
        print(
            "No new hinge comparison written: these groups did not supply "
            "completed results for this run: " + ", ".join(missing_groups)
        )
        return None

    selected_paths = [paths[group].resolve() for group in groups]
    if len(selected_paths) != len(set(selected_paths)):
        raise ValueError("Different groups must have different result workbooks")

    output_file = (
        Path(config["OUTPUT_FOLDER"])
        / "Hinge_Model"
        / "hinge_common_significant.xlsx"
    )
    if output_file.resolve() in selected_paths:
        raise ValueError("The summary cannot overwrite an input result workbook")

    group_tables = []
    group_info = []
    common = None

    for number, group in enumerate(groups, start=1):
        results = _read_hinge_results(paths[group])
        significant = results.loc[
            results["Slope_2_Change_PValue"] < threshold
        ]
        sheet_name = f"Group_{number}_Results"
        group_tables.append((sheet_name, results))
        group_info.append({
            "Group": group,
            "Results_Sheet": sheet_name,
            "Result_File": str(paths[group]),
            "Proteins_Tested": len(results),
            "Significant_Proteins": len(significant),
        })

        # Keep both groups' model details and label columns with actual names.
        labelled = significant.rename(columns={
            column: f"{column}__{group}"
            for column in significant.columns if column != "Protein"
        })
        if common is None:
            common = labelled
        else:
            overlapping_columns = (
                set(common.columns) & set(labelled.columns)
            ) - {"Protein"}
            if overlapping_columns:
                raise ValueError("Group labels produced ambiguous result columns")
            common = common.merge(
                labelled, on="Protein", how="inner", validate="one_to_one"
            )

    common = common.sort_values("Protein").reset_index(drop=True)
    settings = pd.DataFrame([
        {"Setting": "INPUT_FILE_PATH", "Value": str(config.get("INPUT_FILE_PATH", ""))},
        {"Setting": "SHEET_NAME", "Value": str(config.get("SHEET_NAME", ""))},
        {"Setting": "OUTPUT_FOLDER", "Value": str(config["OUTPUT_FOLDER"])},
        {"Setting": "P_THRESHOLD", "Value": threshold},
        {"Setting": "P-value", "Value": "Selected model: Slope_2_Change_PValue (nominal)"},
        {"Setting": "Overlap rule", "Value": "p < P_THRESHOLD in every selected group"},
        {"Setting": "Common proteins", "Value": len(common)},
    ])

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(
        output_file,
        engine="xlsxwriter",
        engine_kwargs={"options": {"strings_to_formulas": False, "strings_to_urls": False}},
    ) as writer:
        common.to_excel(writer, sheet_name="Common_Significant", index=False)
        pd.DataFrame(group_info).to_excel(
            writer, sheet_name="Groups", index=False
        )
        settings.to_excel(writer, sheet_name="Comparison_Info", index=False)
        for sheet_name, results in group_tables:
            results.to_excel(writer, sheet_name=sheet_name, index=False)

    print(
        f"Hinge comparison: {len(common)} proteins with nominal p < {threshold} "
        f"in all {len(groups)} selected groups."
    )
    print(f"Saved: {output_file}")
    return output_file


if __name__ == "__main__":
    raise SystemExit(
        "Run main.py after adding the comparison call at the end of its hinge "
        "loop. This module needs the live configuration and the result paths "
        "returned by that run."
    )

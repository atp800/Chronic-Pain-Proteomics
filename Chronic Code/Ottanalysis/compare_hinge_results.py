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
    required = {
        "Protein",
        "Model_Type",
        "Hinge_Search_PValue",
        "Hinge_FDR",
        "Slope_2_Change_PValue_Nominal"
    }

    missing = required - set(results.columns)
    if missing:
        raise ValueError(
            f"{path}: missing hinge result columns {sorted(missing)}; "
            "rerun hinge analysis."
        )

    results["Protein"] = results["Protein"].str.strip()
    if results["Protein"].fillna("").eq("").any():
        raise ValueError(f"{path}: blank protein identifiers found")
    if results["Protein"].duplicated().any():
        raise ValueError(f"{path}: duplicate protein identifiers found")

    for column in ("Hinge_Search_PValue", "Hinge_FDR", "Slope_2_Change_PValue_Nominal"):
        values = pd.to_numeric(results[column], errors="raise")
        if not values.dropna().between(0, 1).all():
            raise ValueError(f"{path}: {column} must be between 0 and 1")
        results[column] = values
    if (results.Hinge_FDR.notna() & (results.Hinge_Search_PValue.isna()
                                    | results.Hinge_FDR.lt(results.Hinge_Search_PValue))).any():
        raise ValueError(f"{path}: inconsistent search p-values and FDR values")
    fallback = ~results["Model_Type"].eq("Mixed-Effects")

    results.loc[
        fallback,
        [
            "Hinge_Search_PValue",
            "Hinge_FDR",
            "Slope_2_Change_PValue_Nominal"
        ]
    ] = float("nan")

    return results


def _merge_group_overlap(common, significant, group):
    """Intersect a significant protein table with previous groups."""

    labelled = significant.rename(columns={
        column: f"{column}__{group}"
        for column in significant.columns
        if column != "Protein"
    })

    if common is None:
        return labelled

    overlapping_columns = (
        set(common.columns) & set(labelled.columns)
    ) - {"Protein"}

    if overlapping_columns:
        raise ValueError(
            "Group labels produced ambiguous result columns"
        )

    return common.merge(
        labelled,
        on="Protein",
        how="inner",
        validate="one_to_one"
    )


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

    # Existing FDR-significant overlap.
    common = None

    # New nominally significant overlap.
    common_nominal = None

    for number, group in enumerate(groups, start=1):

        results = _read_hinge_results(paths[group])

        # Exclude descriptive OLS fallback models from significance.
        mixed = results["Model_Type"].eq("Mixed-Effects")

        # Original FDR-significant proteins.
        significant_fdr = results.loc[
            mixed & results["Hinge_FDR"].le(threshold)
        ]

        # New nominally significant proteins.
        significant_nominal = results.loc[
            mixed
            & results["Slope_2_Change_PValue_Nominal"].lt(threshold)
        ]

        sheet_name = f"Group_{number}_Results"
        group_tables.append((sheet_name, results))

        group_info.append({
            "Group": group,
            "Results_Sheet": sheet_name,
            "Result_File": str(paths[group]),
            "Proteins_Tested": int(
                results["Hinge_Search_PValue"].notna().sum()
            ),
            "Descriptive_Fallbacks": int(
                results["Model_Type"].eq("Fixed-Effects").sum()
            ),
            # Keep the original field FDR-based for compatibility.
            "Significant_Proteins": len(significant_fdr),
            "Nominal_Significant_Proteins": len(significant_nominal),
        })

        # Intersect FDR-significant proteins across groups.
        common = _merge_group_overlap(
            common,
            significant_fdr,
            group
        )

        # Intersect nominally significant proteins across groups.
        common_nominal = _merge_group_overlap(
            common_nominal,
            significant_nominal,
            group
        )

    # Sort both overlap tables by protein identifier.
    common = common.sort_values(
        "Protein"
    ).reset_index(drop=True)

    common_nominal = common_nominal.sort_values(
        "Protein"
    ).reset_index(drop=True)


    settings = pd.DataFrame([
        {"Setting": "INPUT_FILE_PATH", "Value": str(config.get("INPUT_FILE_PATH", ""))},
        {"Setting": "SHEET_NAME", "Value": str(config.get("SHEET_NAME", ""))},
        {"Setting": "OUTPUT_FOLDER", "Value": str(config["OUTPUT_FOLDER"])},
        {"Setting": "P_THRESHOLD", "Value": threshold},
        {"Setting": "Common_Significant", "Value": "Search-adjusted BH FDR <= P_THRESHOLD in every group"},
        {"Setting": "Common_Nominal", "Value": "AIC-selected hinge nominal slope-change p < P_THRESHOLD in every group; unadjusted for selection/multiple tests"},
        {"Setting": "OLS fallbacks", "Value": "Descriptive only; excluded from both overlaps"},
        {"Setting": "Report ranking mode", "Value": (
                "Nominal"
                if config.get("HINGE_RANK_BY_NOMINAL", False)
                else "FDR"
            )},
        {"Setting": "Interpretation", "Value": "An overlap of groupwise hits is not a direct test of group differences; the nominal overlap is exploratory."},
        {"Setting": "Common proteins (FDR)", "Value": len(common)},
        {"Setting": "Common proteins (nominal)", "Value": len(common_nominal)},
    ])

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(
        output_file,
        engine="xlsxwriter",
        engine_kwargs={"options": {"strings_to_formulas": False, "strings_to_urls": False}},
    ) as writer:
        common.to_excel(
            writer,
            sheet_name="Common_Significant",
            index=False
        )

        common_nominal.to_excel(
            writer,
            sheet_name="Common_Nominal",
            index=False
        )
        pd.DataFrame(group_info).to_excel(
            writer, sheet_name="Groups", index=False
        )
        settings.to_excel(writer, sheet_name="Comparison_Info", index=False)
        for sheet_name, results in group_tables:
            results.to_excel(writer, sheet_name=sheet_name, index=False)

    print(
        f"Hinge comparison: {len(common)} proteins with search-adjusted BH FDR <= {threshold} "
        f"in all {len(groups)} selected groups."
    )
    print(
        f"Hinge comparison: {len(common_nominal)} proteins "
        f"with nominal slope-change p < {threshold} "
        f"in all {len(groups)} selected groups (exploratory)."
    )

    print(f"Saved: {output_file}")
    return output_file


if __name__ == "__main__":
    raise SystemExit(
        "Run main.py after adding the comparison call at the end of its hinge "
        "loop. This module needs the live configuration and the result paths "
        "returned by that run."
    )

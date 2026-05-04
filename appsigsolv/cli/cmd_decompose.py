"""Command logic for 'decompose'."""
import pandas as pd
from pathlib import Path

from appsigsolv.io.data_manager import load_and_preprocess, save_json_config, save_csv
from appsigsolv.core.dia import detect_jumps, auto_detect_periods, prescreen_periods, run_omt_sigma_scan
from appsigsolv.core.modeling import extract_components
from appsigsolv.utils.visualization import save_plot, save_report

def run_decompose(args):
    csv_path = Path(args.input_csv)
    stem = csv_path.stem

    if args.output_dir:
        out_root = Path(args.output_dir) / stem
    else:
        out_root = csv_path.parent / stem
    out_root.mkdir(parents=True, exist_ok=True)

    if args.component.lower() == "all":
        df_peek = pd.read_csv(csv_path, nrows=1)
        detected_date_col = args.date_col
        if not detected_date_col:
            for c in ['gpsdate', 'datetime', 'date', 'Date', 'time', 'Time']:
                if c in df_peek.columns:
                    detected_date_col = c
                    break
            if not detected_date_col:
                detected_date_col = df_peek.columns[0]
                
        components_to_process = [c for c in df_peek.columns if c != detected_date_col]
        print(f"  [batch] Detected {len(components_to_process)} components to process.")
    else:
        components_to_process = [c.strip() for c in args.component.split(",")]

    candidate_periods = [float(p) for p in args.periods.split(",")] if args.periods else []
    extra_jumps = [d.strip() for d in args.jumps.split(",") if d.strip()]
    extra_polylines = [d.strip().replace("-", "") for d in args.polylines.split(",") if d.strip()]
    
    extra_logs = {}
    if args.logs:
        for item in args.logs.split(","):
            if ":" in item:
                dt_str, tau_str = item.split(":", 1)
                dt_key = dt_str.strip().replace("-", "")
                try:
                    tau = float(tau_str)
                    extra_logs[dt_key] = [tau]
                except ValueError:
                    print(f"  [args] Warning: could not parse tau in '{item}'")

    for comp in components_to_process:
        print(f"\n{'='*60}\nProcessing: {comp}\n{'='*60}")

        # Check if already processed
        json_path = out_root / f"{stem}_model_{comp}.json"
        if json_path.exists():
            print(f"  [info] Component '{comp}' already processed. Skipping.")
            continue

        try:
            series, detected_date_col = load_and_preprocess(str(csv_path), comp, args.date_col, args.unit, irregular=args.irregular)
        except KeyError:
            print(f"  [error] Column '{comp}' not found. Skipping.")
            continue
            
        print(f"  Series loaded: {len(series)} points")

        if len(series) == 0:
            print(f"  [error] Component '{comp}' has no valid data after preprocessing. Skipping.")
            continue

        jump_dates = detect_jumps(series, extra_jumps)
        
        final_periods = list(candidate_periods)
        if final_periods:
            print(f"  [periods] Forcing user-provided periods: {final_periods}")
            
        if args.auto_periods > 0:
            auto_p = auto_detect_periods(series, max_periods=args.auto_periods)
            screened_auto = prescreen_periods(series, auto_p)
            for p in screened_auto:
                if not any(abs(p - existing) < 0.05 for existing in final_periods):
                    final_periods.append(p)
                    
        final_periods.sort()

        candidate_degrees = [0, 1, 2] if args.poly_deg == -1 else [args.poly_deg]
        best_overall_model = None
        best_overall_scan_table = []

        for deg in candidate_degrees:
            if args.poly_deg == -1:
                print(f"\n  [auto-deg] Testing polynomial degree: {deg}")

            best_model_deg, scan_table_deg = run_omt_sigma_scan(
                series, jump_dates, final_periods, extra_polylines, extra_logs, deg,
                args.sigma_min, args.sigma_max, args.sigma_step,
                args.alpha, args.max_iter,
                no_relax=args.no_relax, cores=args.cores
            )
            
            if best_model_deg is not None:
                if best_overall_model is None:
                    best_overall_model = best_model_deg
                    best_overall_scan_table = scan_table_deg
                else:
                    curr_stats = best_overall_model["_omt_stats"]
                    new_stats = best_model_deg["_omt_stats"]
                    if (new_stats["n_param"] < curr_stats["n_param"]) or \
                       (new_stats["n_param"] == curr_stats["n_param"] and new_stats["p_value"] > curr_stats["p_value"]):
                        best_overall_model = best_model_deg
                        best_overall_scan_table = scan_table_deg

        if best_overall_model is None:
            print(f"  WARNING: No accepted model found for {comp} (tried degrees {candidate_degrees}). Skipping output.")
            continue

        best_model = best_overall_model
        scan_table = best_overall_scan_table
        
        final_deg = best_model.get("polynomial", 0)
        print(f"  [auto-deg] Selected polynomial degree: {final_deg}")

        extracted = extract_components(series, best_model, comp)
        print(f"  [extract] Components: {list(extracted.keys())}")
        
        save_json_config(best_model, comp, out_root, stem)
        save_csv(series, extracted, comp, out_root, stem)
        if not args.no_plot:
            save_plot(series, extracted, best_model, comp, args.unit, out_root, stem)
        save_report(series, extracted, best_model, scan_table, comp, args.unit, out_root, stem)

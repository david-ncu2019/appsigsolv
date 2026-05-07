"""Command logic for 'reconstruct'."""
import os
import calendar
import pandas as pd
from appsigsolv.io.data_manager import load_and_clean_data_for_reconstruct, load_json_config
from appsigsolv.core.modeling import estimate_time_func, get_design_matrix4time_func


def _build_custom_dates(start, end, day_list):
    """Generate Timestamps for given days-of-month between start and end, clamped to month length."""
    result = []
    current = start.replace(day=1)
    while current <= end:
        last_day = calendar.monthrange(current.year, current.month)[1]
        for d in sorted(day_list):
            actual_day = min(d, last_day)
            ts = pd.Timestamp(current.year, current.month, actual_day)
            if start <= ts <= end:
                result.append(ts)
        # Advance to next month
        if current.month == 12:
            current = current.replace(year=current.year + 1, month=1)
        else:
            current = current.replace(month=current.month + 1)
    return pd.DatetimeIndex(result)

def run_reconstruct(args):
    try:
        df_in, date_col, target_col, dates, raw_disp = load_and_clean_data_for_reconstruct(
            args.input_file, args.date_col, args.target_col, args.unit
        )
    except Exception as e:
        print(f"Error: {e}")
        return

    model = {
        'polynomial': 1,
        'periodic': [],
        'stepDate': [],
        'polyline': [],
        'exp': {},
        'log': {},
        'exp_trend': None,
    }
    
    if args.json_file:
        json_config = load_json_config(args.json_file)
        model.update(json_config)
        
    if args.poly is not None: model['polynomial'] = args.poly
    if args.period is not None: model['periodic'] = args.period
    if args.stepDate is not None: model['stepDate'] = args.stepDate
    if args.polyline is not None: model['polyline'] = args.polyline
    
    if args.exp:
        for onset, tau in args.exp:
            onset_key = onset.replace('-', '')
            if onset_key not in model['exp']: model['exp'][onset_key] = []
            model['exp'][onset_key].append(float(tau))
            
    if args.log:
        for onset, tau in args.log:
            onset_key = onset.replace('-', '')
            if onset_key not in model['log']: model['log'][onset_key] = []
            model['log'][onset_key].append(float(tau))

    exp_trend_arg = getattr(args, "exp_trend", None)
    if exp_trend_arg is not None:
        model['exp_trend'] = float(exp_trend_arg)

    print(f"Fitting model: {model}")
    try:
        ref_date = args.ref_date if args.ref_date else None
        disp_ref = raw_disp - raw_disp[0]
        
        G, m_est, e2, d_hat = estimate_time_func(model, dates, disp_ref)
        
        if args.unit == 'mm':
            modeled_out = (d_hat + raw_disp[0]) * 1000.0
        else:
            modeled_out = d_hat + raw_disp[0]
            
    except Exception as e:
        print(f"Error fitting model: {e}")
        return

    if args.sampling_rate == 'daily':
        print("Generating daily reconstruction...")
        out_dates = pd.date_range(start=min(dates), end=max(dates), freq='D')

    elif args.sampling_rate == 'custom':
        if not args.custom_dates:
            print(
                "Error: --sampling-rate custom requires --custom-dates.\n"
                "  Example: --custom-dates 1,6,11,16,21,26\n"
                "  Provide comma-separated day numbers (1–31). Days exceeding the month length\n"
                "  are clamped to the last day of that month (e.g., 31 → Feb 28/29)."
            )
            return

        try:
            day_list = [int(d.strip()) for d in args.custom_dates.split(',')]
        except ValueError:
            print(
                "Error: --custom-dates must be comma-separated integers.\n"
                "  Example: --custom-dates 1,6,11,16,21,26"
            )
            return

        invalid = [d for d in day_list if not (1 <= d <= 31)]
        if invalid:
            print(
                f"Error: Day values out of range (1–31): {invalid}\n"
                "  Example: --custom-dates 1,6,11,16,21,26"
            )
            return

        print(f"Generating custom reconstruction at days: {day_list} of each month...")
        out_dates = _build_custom_dates(min(dates), max(dates), day_list)

    else:
        out_dates = None

    if out_dates is not None:
        date_list = [d.to_pydatetime() for d in out_dates]
        G_out = get_design_matrix4time_func(date_list, model, ref_date=dates[0])
        d_out = G_out @ m_est

        if args.unit == 'mm':
            d_out_final = (d_out + raw_disp[0]) * 1000.0
        else:
            d_out_final = d_out + raw_disp[0]

        df_out = pd.DataFrame({
            date_col: out_dates,
            'reconstructed': d_out_final
        })
    else:
        df_out = df_in.copy()
        df_out['modeled'] = modeled_out

    if not args.outfile:
        base, ext = os.path.splitext(args.input_file)
        suffix = f"_{args.sampling_rate}" if args.sampling_rate else "_modeled"
        args.outfile = f"{base}{suffix}{ext}"
    
    print(f"Saving results to {args.outfile}...")
    if args.outfile.endswith('.csv'):
        df_out.to_csv(args.outfile, index=False)
    else:
        df_out.to_excel(args.outfile, index=False)

    print("Done.")

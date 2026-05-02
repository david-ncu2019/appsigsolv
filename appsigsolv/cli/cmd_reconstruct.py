"""Command logic for 'reconstruct'."""
import os
import pandas as pd
from appsigsolv.io.data_manager import load_and_clean_data_for_reconstruct, load_json_config
from appsigsolv.core.modeling import estimate_time_func, get_design_matrix4time_func

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
        'log': {}
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

    if args.daily:
        print("Generating daily reconstruction...")
        start_date = min(dates)
        end_date = max(dates)
        daily_dates = pd.date_range(start=start_date, end=end_date, freq='D')
        daily_list = [d.to_pydatetime() for d in daily_dates]
        
        G_daily = get_design_matrix4time_func(daily_list, model, ref_date=dates[0])
        d_daily = G_daily @ m_est
        
        if args.unit == 'mm':
            d_daily_out = (d_daily + raw_disp[0]) * 1000.0
        else:
            d_daily_out = d_daily + raw_disp[0]
            
        df_out = pd.DataFrame({
            date_col: daily_dates,
            'reconstructed': d_daily_out
        })
    else:
        df_out = df_in.copy()
        df_out['modeled'] = modeled_out

    if not args.outfile:
        base, ext = os.path.splitext(args.input_file)
        suffix = "_daily" if args.daily else "_modeled"
        args.outfile = f"{base}{suffix}{ext}"
    
    print(f"Saving results to {args.outfile}...")
    if args.outfile.endswith('.csv'):
        df_out.to_csv(args.outfile, index=False)
    else:
        df_out.to_excel(args.outfile, index=False)

    print("Done.")

"""Argparse configuration and dispatcher."""
import argparse
from appsigsolv.cli.cmd_decompose import run_decompose
from appsigsolv.cli.cmd_reconstruct import run_reconstruct

def create_parser():
    parser = argparse.ArgumentParser(
        description="appsigsolv - Applied Geology's Signal Solver for Timeseries Deformation",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")
    subparsers.required = True

    # ---------------- DECOMPOSE ----------------
    p_dec = subparsers.add_parser("decompose", help="GPS signal decomposition using OMT-based parametric modeling")
    p_dec.add_argument("input_csv", help="Path to input CSV file")
    p_dec.add_argument("--component", default="dU", help="Component(s) to process: e.g., dU | all (default: dU)")
    p_dec.add_argument("--date-col", default="", help="Specify the date column name")
    p_dec.add_argument("--unit", choices=['m', 'mm'], default="mm", help="Input data unit (default: mm)")
    p_dec.add_argument("--jumps", default="", help="Extra jump dates YYYY-MM-DD,YYYY-MM-DD")
    p_dec.add_argument("--polylines", default="", help="Extra polyline break dates YYYY-MM-DD,YYYY-MM-DD")
    p_dec.add_argument("--logs", default="", help="Extra log relaxation dates YYYY-MM-DD:tau,YYYY-MM-DD:tau")
    p_dec.add_argument("--poly-deg", type=int, default=1, help="Polynomial degree for trend: 0=offset, 1=linear, 2=accel (default: 1)")
    p_dec.add_argument("--periods", default="0.25,0.5,1.0,2.0", help="Candidate periods in years (default: 0.25,0.5,1.0,2.0)")
    p_dec.add_argument("--auto-periods", type=int, default=5, help="Auto-detect up to N additional dominant periods (default: 5)")
    p_dec.add_argument("--sigma-min", type=float, default=2.0, help="Min sigma (default 2.0)")
    p_dec.add_argument("--sigma-max", type=float, default=15.0, help="Max sigma (default 15.0)")
    p_dec.add_argument("--sigma-step", type=float, default=0.5, help="Sigma step (default 0.5)")
    p_dec.add_argument("--alpha", type=float, default=0.05, help="OMT significance level (default 0.05)")
    p_dec.add_argument("--max-iter", type=int, default=5, help="Max DIA iterations per sigma (default 5)")
    p_dec.add_argument("--irregular", action="store_true", help="Do NOT resample to a daily grid")
    p_dec.add_argument("--no-plot", action="store_true", help="Skip PNG generation")
    p_dec.add_argument("--no-relax", action="store_true", help="Skip exponential/logarithmic relaxation testing")
    p_dec.add_argument("--output-dir", default="", help="Parent dir for output folder")
    p_dec.add_argument("--cores", type=int, default=1, help="Number of CPU cores for parallel processing (default: 1)")
    p_dec.add_argument(
        "--exp-trend", dest="exp_trend", default=None, metavar="B_OR_AUTO",
        help="Exponential trend component exp(-b*t)-1. Pass 'auto' to auto-detect the best decay rate, "
             "or a numeric b value in 1/days (e.g. 0.001 ≈ 2.7-yr time constant). Omit to disable."
    )
    
    # ---------------- RECONSTRUCT ----------------
    p_rec = subparsers.add_parser("reconstruct", help="Reconstruct 1D timeseries signal using parametric deformation modeling")
    p_rec.add_argument("input_file", help="Input CSV or Excel file")
    p_rec.add_argument("--json", dest="json_file", help="Path to JSON configuration file")
    p_rec.add_argument("--date-col", help="Date column name")
    p_rec.add_argument("--target-col", help="Target displacement column name")
    p_rec.add_argument("--unit", choices=['m', 'mm'], default='mm', help="Input unit (default: mm)")
    p_rec.add_argument("--poly", type=int, help="Polynomial degree (0=offset, 1=velocity, 2=accel)")
    p_rec.add_argument("--period", type=float, nargs="*", help="Periodic components in years")
    p_rec.add_argument("--step", dest="stepDate", nargs="*", help="Step dates in YYYYMMDD")
    p_rec.add_argument("--polyline", nargs="*", help="Polyline break dates in YYYYMMDD")
    p_rec.add_argument("--exp", nargs=2, action="append", metavar=("DATE", "TAU"), help="Exponential: onset (YYYYMMDD) and tau (days)")
    p_rec.add_argument("--log", nargs=2, action="append", metavar=("DATE", "TAU"), help="Logarithmic: onset (YYYYMMDD) and tau (days)")
    p_rec.add_argument(
        "--exp-trend", dest="exp_trend", default=None, type=float, metavar="B",
        help="Exponential trend b value in 1/days for exp(-b*t)-1 component (e.g. 0.001). "
             "Loaded automatically from --json if present."
    )
    p_rec.add_argument("-o", "--output", dest="outfile", help="Output file name")
    p_rec.add_argument(
        "--sampling-rate",
        dest="sampling_rate",
        choices=["daily", "custom"],
        help="Output sampling: 'daily' for every day, 'custom' for specific days of each month"
    )
    p_rec.add_argument(
        "--custom-dates",
        dest="custom_dates",
        help="Comma-separated days of month for --sampling-rate custom (e.g. 1,6,11,16,21,26)"
    )
    p_rec.add_argument("--ref-date", help="Reference date for model (YYYYMMDD)")
    
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    
    if args.command == "decompose":
        run_decompose(args)
    elif args.command == "reconstruct":
        run_reconstruct(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

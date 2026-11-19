#!/usr/bin/env bash
# Stubs only the processor invocation: no real GCOV data is required.
set -euo pipefail
basedir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/input"
touch "$tmp/input/NISAR_L2_PR_GCOV_001_001_A_001_2005_TEST_A_20261119T000000_20261119T000100_P05036_P_F_A_001.h5"
cat > "$tmp/fake_processor.py" <<'PY'
import pathlib, sys
args = sys.argv[1:]
for option in ('--input_dir','--output_dir','--ease2_res','--backscatter','--filter'):
    assert option in args, (option,args)
assert args[args.index('--ease2_res')+1] == '50'
assert args[args.index('--backscatter')+1] == 'sigma0'
output = pathlib.Path(args[args.index('--output_dir')+1])
assert output.is_dir(), output
assert list(pathlib.Path(args[args.index('--input_dir')+1]).glob('*.h5'))
(output/'test_success.txt').write_text('ok')
PY
(
 cd "$tmp"
 MAAP_PYTHON_EXEC="$(command -v python3)" \
 MAAP_PROCESSOR_PATH="$tmp/fake_processor.py" \
 bash "$basedir/run_maap.sh" \
   --input_file "https://example.test/$(basename "$tmp/input"/*.h5)" \
   --ease2_res 50 --backscatter sigma0 --mode all
)
for f in nofilt enhanced_lee gamma_map; do
  test -s "$tmp/output/$f/test_success.txt"
done
echo 'PASS: MAAP foreground wrapper resolves staged URLs, runs 3 modes, writes output/'

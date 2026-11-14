"""Amendment 1 runner for WP3 (see AMEND-WP3-1.md). Imports the sealed wp3.py unchanged and changes ONLY the two
bulk refinement-stopping tolerances before calling its main(). Usage:
  memcap 8 py -3.11 wp3_a1.py --main result_wp3.json   (exit 3 = rerun the same command)"""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wp3
wp3.BULK_TOLF = 3e-2    # was 1e-3; 3e-2 is the value the sealed file's own controls use (line 2738)
wp3.RCAP_FRAC = 1.0     # was 0.8; bulk Ricci boxes may stop at the neck's leading value L (a lower bound for R_out)
if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != '--main': sys.exit('usage: wp3_a1.py --main OUT.json')
    print('AMEND-WP3-1: BULK_TOLF = %g, RCAP_FRAC = %g' % (wp3.BULK_TOLF, wp3.RCAP_FRAC), flush=True)
    wp3.main(sys.argv[2], 530.0)

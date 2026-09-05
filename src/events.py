"""Small FITS event-file utilities used by the pipeline scripts (numpy + astropy)."""

import hashlib

import numpy as np
from astropy.io import fits
from astropy.table import Table


def read_events(path):
    """Open an event file; return (hdul, events_hdu, gti) with gti an (N, 2) float array."""
    hdul = fits.open(path)
    events = hdul['EVENTS']
    gti_hdu = hdul['GTI']
    gti = np.column_stack([gti_hdu.data['START'].astype(float), gti_hdu.data['STOP'].astype(float)])
    return hdul, events, gti


def check_time_header(hdr):
    """Assert that the file is barycentered (TDB, SOLARSYSTEM, TIMEZERO 0); return the keywords."""
    timezero = float(hdr.get('TIMEZERO', 0.0))
    assert timezero == 0.0, f'TIMEZERO = {timezero}, expected 0'
    assert str(hdr['TIMESYS']).startswith('TDB'), f"TIMESYS = {hdr['TIMESYS']}, expected TDB"
    assert str(hdr['TIMEREF']).strip() == 'SOLARSYSTEM', f"TIMEREF = {hdr['TIMEREF']}, expected SOLARSYSTEM"
    return {'TIMESYS': str(hdr['TIMESYS']).strip(), 'TIMEREF': str(hdr['TIMEREF']).strip(), 'TIMEZERO': timezero,
            'TELESCOP': str(hdr.get('TELESCOP', '')).strip(), 'INSTRUME': str(hdr.get('INSTRUME', '')).strip()}


def md5sum(path):
    """MD5 hex digest of a file."""
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def energy_filter_events(in_path, out_path, pi_lo, pi_hi, history):
    """Keep events with pi_lo <= PI <= pi_hi (inclusive); copy every other HDU. Returns the count kept."""
    with fits.open(in_path) as hdul:
        data = hdul['EVENTS'].data
        pi = np.asarray(data['PI'])
        mask = (pi >= pi_lo) & (pi <= pi_hi)
        hdul['EVENTS'].data = data[mask]
        for line in history:
            hdul['EVENTS'].header.add_history(line)
        hdul.writeto(out_path, overwrite=True, checksum=True)
    return int(mask.sum())


def apply_gti(in_path, out_path, gti, history):
    """Keep events inside the GTI table (start <= t <= stop), replace the GTI HDU, update EXPOSURE/ONTIME."""
    with fits.open(in_path) as hdul:
        data = hdul['EVENTS'].data
        t = np.asarray(data['TIME'], dtype=np.float64)

        # searchsorted on the sorted GTI starts; an event is inside when it lies before that GTI's stop
        start = gti[:, 0]
        stop = gti[:, 1]
        idx = np.searchsorted(start, t, side='right') - 1
        inside = (idx >= 0) & (t <= stop[np.clip(idx, 0, len(stop) - 1)])
        hdul['EVENTS'].data = data[inside]

        new_gti = fits.BinTableHDU.from_columns([
            fits.Column(name='START', format='D', unit='s', array=start),
            fits.Column(name='STOP', format='D', unit='s', array=stop)], header=hdul['GTI'].header)
        hdul[hdul.index_of('GTI')] = new_gti

        exposure = float(np.sum(stop - start))
        for name in ['EVENTS', 'GTI']:
            hdul[name].header['EXPOSURE'] = exposure
            hdul[name].header['ONTIME'] = exposure
        for line in history:
            hdul['EVENTS'].header.add_history(line)
        hdul.writeto(out_path, overwrite=True, checksum=True)
    return int(inside.sum())


def read_gti_table(path):
    """Read a GTI table saved as .ecsv; return (gti array (N, 2), metadata dict)."""
    table = Table.read(path, format='ascii.ecsv')
    gti = np.column_stack([np.asarray(table['START'], dtype=np.float64), np.asarray(table['STOP'], dtype=np.float64)])
    return gti, dict(table.meta)


def parse_tim(path):
    """Parse a tempo2-format .tim file into a list of dicts (the MJD string is kept verbatim)."""
    toas = []
    with open(path) as f:
        for raw in f:
            line = raw.rstrip('\n')
            if line.strip() == '' or line.startswith('FORMAT') or line.startswith('C ') or line.startswith('#'):
                continue
            parts = line.split()
            flags = {parts[k].lstrip('-'): parts[k + 1] for k in range(5, len(parts) - 1, 2)}
            toas.append({'mjd': parts[2], 'err_us': float(parts[3]), 'flags': flags, 'line': line})
    return toas


def write_tim(path, toas, comments=None, extra_flags=None):
    """Write TOA dicts (from parse_tim) as a tempo2 .tim file with optional 'C ' comment lines."""
    with open(path, 'w') as f:
        f.write('FORMAT 1\n')
        for c in (comments or []):
            f.write(f'C {c}\n')
        for toa in toas:
            line = toa['line']
            if extra_flags:
                line = line + ' ' + ' '.join(f'-{k} {v}' for k, v in extra_flags.items())
            f.write(line + '\n')

"""
Access Lawrence Berkeley National Laboratory's X-ray reflectivity database for index of refraction coefficents as a function of energy
"""

# Standard library modules.
import urllib

# Third party modules.
import requests
import bs4
import numpy as np
import time

# Local modules.
import pdb

# Globals and constants variables.

class _script_call(object):

    URLBASE = 'http://henke.lbl.gov'
    MAX_STEP = 500

    def _process(self, script, data):
        response = self._post(script, data)
        self._check_errors(response)
        s = self._retrieve_data(response)
        return self._parse_data(s)

    def _post(self, script, data):
        url = urllib.request.urljoin(self.URLBASE, '/cgi-bin/' + script)
        #time.sleep(10)
        return requests.post(url, data)

    def _check_errors(self, response):
        if not response.status_code == requests.codes.ok: #@UndefinedVariable
            raise RuntimeError('Could not connect to server')

        #if not response.text.startswith('<Head>'):
        #    raise RuntimeError(response.text)

    def _retrieve_data(self, response):
        #soup = bs4.BeautifulSoup(response.text, 'html.parser')
        #urldata = soup.body.img.attrs['src'].replace('.gif','.dat')
        urldata = response.text[response.text.find('HREF') + 6:response.text.find('here') - 2]
        #pdb.set_trace()
        url = urllib.request.urljoin(self.URLBASE, urldata)
        with urllib.request.urlopen(url) as fp:
            return fp.read()

    def _parse_data(self, s):
        return [list(map(float, line.split())) for line in s.splitlines()[2:]]

    def _iter_range(self, x0, x1, step):
        dx = (x1 - x0) / step
        xs = np.arange(x0, x1, dx)
        x0s = xs[::self.MAX_STEP]
        x1s = x0s + dx * min(step, self.MAX_STEP)
        x1s[-1] = x1
        steps = np.array((x1s - x0s) / dx, dtype=int)
        yield from zip(x0s, x1s, steps)

    #def calculate_energy_scan(self, e0_eV, e1_eV, step):
    #    """
    #    Returns three columns:
        
    #      * Photon Energy (eV)
    #      * Reflectivity
    #      * Transmission into substrate
    #    """
    #    raise NotImplementedError

    #def calculate_wavelength_scan(self, lambda0_nm, lambda1_nm, step):
    #    """
    #    Returns three columns:
        
    #      * Wavelength (nm)
    #      * Reflectivity
    #      * Transmission into substrate
    #    """
    #    raise NotImplementedError

    #def calculate_angle_scan(self, theta0_deg, theta1_deg, step, energy_eV):
    #    """
    #    Returns three columns:
        
    #      * Angle (deg)
    #      * Reflectivity
    #      * Transmission into substrate
    #    """
    #    raise NotImplementedError

class index(_script_call):

    def __init__(self):
        self.chemical_formula = 'Si3N4'
        self.density_g_cm3 = -1

    def _create_post_data(self):
        data = {}
        data['Formula'] = self.chemical_formula
        data['Density'] = self.density_g_cm3
        data['Output'] = 'Text File'
        return data

    def calculate_energy_scan(self, e0_eV, e1_eV, step = 500):
        data = self._create_post_data()
        data['Scan'] = 'Energy'
        data['Min'] = e0_eV
        data['Max'] = e1_eV
        data['Npts'] = step
        return self._process('getdb.pl', data)
import pyvo as vo
import csv
import pathlib
import tarfile

from astropy.utils.data import download_file
import matplotlib.pyplot as plt
import shutil
from astroquery.eso import Eso
import pandas as pd
import requests
from astropy.table import Table
from astropy.io import fits
import os
import numpy as np
from scipy.signal import find_peaks
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
##import dustmaps
##from dustmaps.config import config
##import dustmaps.bayestar

def downloading_data():
    from astroquery.gaia import Gaia
    cone_radius = 0.1

    gaia_service = vo.dal.TAPService(
        "https://gea.esac.esa.int/tap-server/tap"
    )

    adql_query = """
    select
        bm.source_id,
        gs.ra,
        gs.dec,
        bm.m1,
        bm.m1_lower,
        bm.m1_upper,
        bm.m2,
        bm.m2_lower,
        bm.m2_upper,
        nss.nss_solution_type,
        nss.significance,
        nss.period,
        nss.eccentricity,
        nss.semi_amplitude_primary,
        nss.center_of_mass_velocity,
        gs.ruwe,
        gs.radial_velocity,
        nss.goodness_of_fit
    from gaiadr3.binary_masses as bm
    join gaiadr3.nss_two_body_orbit as nss
        on bm.source_id = nss.source_id
    join gaiadr3.gaia_source as gs
        on bm.source_id = gs.source_id
    where bm.m2_lower > bm.m1_upper
      and bm.m2_lower > 1.5
      and bm.m2_lower < 3
      and nss.nss_solution_type = 'SB1'
      and nss.significance > 20
    """

    query_results = gaia_service.search(adql_query)
    results_table = query_results.to_table()

    print(results_table)
    print(f"Number of candidates: {len(results_table)}")


    with open('Gaia_Data_DR3.csv', 'w', newline='') as f:
        writer = csv.writer(f)

        writer.writerow(results_table.colnames)
        writer.writerows(results_table.iterrows())

    eso = Eso()

    eso_service = vo.dal.TAPService("https://archive.eso.org/tap_obs")

    uves_table = []

    for target in results_table:
        ra = target["ra"]
        dec = target["dec"]
        source_id = target["source_id"]

        adql_query = f"""
        select *
        from ivoa.ObsCore
        where instrument_name = 'UVES'
          and intersects(
              circle('ICRS', {ra}, {dec}, {cone_radius}),
              s_region
          ) = 1
        """

        eso_query = eso_service.search(adql_query).to_table()

        print("UVES", source_id, len(eso_query))

        if len(eso_query) > 0:
            eso_query["GAIA_SOURCE_ID"] = source_id
            uves_table.append(eso_query)

    with open('UVES_data.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(uves_table[0].colnames)

        for table in uves_table:
            writer.writerows(table)


    feros_table = []

    for target in results_table:
        ra = target["ra"]
        dec = target["dec"]
        source_id = target["source_id"]

        adql_query = f"""
        select *
        from ivoa.ObsCore
        where instrument_name = 'FEROS'
          and intersects(
              circle('ICRS', {ra}, {dec}, {cone_radius}),
              s_region
          ) = 1
        """

        eso_query = eso_service.search(adql_query).to_table()

        print("FEROS", source_id, len(eso_query))

        if len(eso_query) > 0:
            eso_query["GAIA_SOURCE_ID"] = source_id
            feros_table.append(eso_query)

    with open('FEROS_data.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(feros_table[0].colnames)

        for table in feros_table:
            writer.writerows(table)


    datatype = "both"

    if datatype == "feros":


        os.makedirs("feros_data", exist_ok=True)

        for target in feros_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"feros_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"feros_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)

    elif datatype == "uves":

        os.makedirs("uves_data", exist_ok=True)

        for target in uves_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"uves_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"uves_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)



    elif datatype == "both":
        os.makedirs("uves_data", exist_ok=True)

        for target in uves_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"uves_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"uves_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)



        os.makedirs("feros_data", exist_ok=True)

        for target in feros_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"feros_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"feros_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)
        



def barycentric_corrections():

    uves_array = os.listdir("uves_data")
    paranal = EarthLocation.of_site("paranal")
    for folder in uves_array:
        [os.remove(f"uves_data/{folder}/{file}") for file in os.listdir(f"uves_data/{folder}") if file.endswith("_bary.fits")]
        folder_data = os.listdir(f"uves_data/{folder}")




        UVES_csv = pd.read_csv("UVES_data.csv")

        
    

        row = UVES_csv[UVES_csv["GAIA_SOURCE_ID"] == int(folder)].iloc[0]
        ra = row["s_ra"]
        dec = row["s_dec"]

        
        coord = SkyCoord(ra * u.deg, dec * u.deg)

        
        for file in folder_data:
            hdul = fits.open(f"uves_data/{folder}/{file}")
            header = hdul[0].header
            
  
            
            if header["SPECSYS"] == "TOPOCENT":

                
                print(f"{file} is topocentric, performing barycentric corrections")

                name = file.replace(".fits", "")
                shutil.copy(f"uves_data/{folder}/{name}.fits", f"uves_data/{folder}/{name}_bary.fits")

                hdul = fits.open(f"uves_data/{folder}/{name}_bary.fits", mode="update")
                data = hdul[1].data

                wavelength = data["WAVE"][0]
                obstime = Time(hdul[1].header["TMID"], format="mjd", scale="utc")
                barycor = coord.radial_velocity_correction(kind="barycentric", obstime=obstime, location=paranal).to(u.km / u.s).value

                
                
                print(f"Barycor: {barycor}")



                data["WAVE"][0] = wavelength * (1 + barycor / 299792.458)
                hdul[0].header["SPECSYS"] = "BARYCENT"

                hdul.close()






                
            elif header["SPECSYS"] == "BARYCENT":
                print(f"{file} is barycentric, skipping")
                shutil.copy(f"{file}.fits", f"{file}_bary.fits")


                
            else:
                print("SYSSPEC not recognised, skipping")
                continue
        


    feros_array = os.listdir("feros_data")

    for folder in feros_array:
        [os.remove(f"feros_data/{folder}/{file}") for file in os.listdir(f"feros_data/{folder}") if file.endswith("_bary.fits")]
        folder_data = os.listdir(f"feros_data/{folder}")

        for file in folder_data:
            file_path = f"feros_data/{folder}/{file}"

            with fits.open(file_path) as hdul:
                header = hdul[0].header
                specsyst = header.get("SPECSYS")

            if specsyst == "TOPOCENT":
                print(f"{file} is topocentric, performing barycentric corrections")

                name = file.replace(".fits", "")
                shutil.copy(f"feros_data/{folder}/{name}.fits", f"feros_data/{folder}/{name}_bary.fits")

                hdul = fits.open(f"feros_data/{folder}/{name}_bary.fits", mode="update")
                data = hdul[1].data

                wavelength = data["WAVE"][0]
                barycor = header["ESO DRS BARYCORR"]

                data["WAVE"][0] = wavelength * (1 + barycor / 299792.458)
                hdul[0].header["SPECSYS"] = "BARYCENT"

                hdul.close()









            elif specsyst == "BARYCENT":
                print(f"{file} is barycentric, copying file")

                name, ext = os.path.splitext(file)
                new_path = f"feros_data/{folder}/{name}_bary{ext}"

                shutil.copy(file_path, new_path)

            else:
                print(f"{file}: SPECSYS not recognised, skipping")
                continue
            




def plot_spectra_old():

    uves_array = os.listdir("uves_data")
    os.makedirs("spectra_UVES_old", exist_ok=True)

    os.makedirs("combined_spectra_old", exist_ok=True)
    os.makedirs("combined_spectra_old/FEROS", exist_ok=True)
    os.makedirs("combined_spectra_old/UVES", exist_ok=True)

    feros_array = os.listdir("feros_data")
    os.makedirs("spectra_FEROS_old", exist_ok=True)

    grid_spectra_uves = []
    grid_spectra_feros = []

    if len(uves_array) == 0:
        print("No UVES data found. Skipping")

    else:


        print(uves_array)
        peak_rows = []
        for folder in uves_array:

            folder_data = [file for file in os.listdir(f"uves_data/{folder}") if not file.endswith("_bary.fits")]
            os.makedirs(f"spectra_UVES_old/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:

                hdul = fits.open(f"uves_data/{folder}/{file}")
                data = hdul[1].data

                wavelength = data["WAVE"][0]
                flux = data["FLUX"][0]


                plt.plot(data["WAVE"][0], data["FLUX"][0], linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)")
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_UVES_old/{folder}/{file}_spectra_UVES_old.png")
                plt.close()
                print(f"(old) UVES: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(data["FLUX"][0])


            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_uves.append((folder, common_wave, average_flux))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)")
            plt.savefig(f"combined_spectra_old/UVES/average_spectrum_{folder}_old.png")
            plt.close()
            print(f"(old) UVES: Combined spectra {folder} saved.")



    if len(feros_array) == 0:
        print("No feros data found. Skipping")

    else:


        print(feros_array)
        for folder in feros_array:

            folder_data = [file for file in os.listdir(f"feros_data/{folder}") if not file.endswith("_bary.fits")]
            os.makedirs(f"spectra_FEROS_old/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:

                hdul = fits.open(f"feros_data/{folder}/{file}")
                data = hdul[1].data



                plt.plot(data["WAVE"][0], data["FLUX"][0], linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux (ADU)")
                plt.ylim(np.nanpercentile(data["FLUX"][0], [1, 99]))
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_FEROS_old/{folder}/{file}_spectra_FEROS_old.png")
                plt.close()
                print(f"(old) FEROS: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(data["FLUX"][0])


            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_feros.append((folder, common_wave, average_flux))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux (ADU)")
            plt.ylim(np.nanpercentile(average_flux, [1, 99]))
            plt.savefig(f"combined_spectra_old/FEROS/average_spectrum_{folder}_old.png")
            plt.close()
            print(f"(old) FEROS: Combined spectra {folder} saved.")





    flux_labels = {
        "UVES": "Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)",
        "FEROS": "Flux (ADU)"
    }

    wave_ranges = {
        "UVES": [("blue", 0, 5760), ("red", 5760, np.inf)],
        "FEROS": [("full", 0, np.inf)]
    }

    for instrument, grid_spectra in [
        ("UVES", grid_spectra_uves),
        ("FEROS", grid_spectra_feros)
    ]:

        for range_name, wave_lo, wave_hi in wave_ranges[instrument]:

            masks = [(wave >= wave_lo) & (wave < wave_hi) & (flux != 0) for _, wave, flux in grid_spectra]

            x_min = min(np.min(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            x_max = max(np.max(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_min = min(np.nanpercentile(flux[mask], 1) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_max = max(np.nanpercentile(flux[mask], 99) for (_, wave, flux), mask in zip(grid_spectra, masks))

            for start in range(0, len(grid_spectra), 8):

                group = grid_spectra[start:start + 8]
                nrows = (len(group) + 1) // 2

                fig, axes = plt.subplots(nrows, 2, figsize=(14, 3 * nrows), sharex=True, sharey=True, squeeze=False)
                axes = axes.flatten()

                for ax, (folder, wavelength, flux) in zip(axes, group):
                    ax.plot(wavelength, flux, linewidth=0.7)
                    ax.set_title(folder, fontsize=9)
                    ax.set_xlim(x_min, x_max)
                    ax.set_ylim(y_min, y_max)
                    ax.grid(True, alpha=0.2)

                for ax in axes[len(group):]:
                    ax.set_visible(False)

                fig.supxlabel("Wavelength (Å)")
                fig.supylabel(flux_labels[instrument])
                fig.tight_layout()

                fig.savefig(
                    f"combined_spectra_old/{instrument}_spectral_grid_{range_name}_{start // 8 + 1}_old.png",
                    dpi=300,
                    bbox_inches="tight"
                )
                plt.close(fig)













def plot_spectra_bary():

    uves_array = os.listdir("uves_data")
    os.makedirs("spectra_UVES_bary", exist_ok=True)

    os.makedirs("combined_spectra_bary", exist_ok=True)
    os.makedirs("combined_spectra_bary/FEROS", exist_ok=True)
    os.makedirs("combined_spectra_bary/UVES", exist_ok=True)

    feros_array = os.listdir("feros_data")
    os.makedirs("spectra_FEROS_bary", exist_ok=True)

    grid_spectra_uves = []
    grid_spectra_feros = []

    if len(uves_array) == 0:
        print("No UVES data found. Skipping")

    else:


        print(uves_array)
        for folder in uves_array:

            folder_data = [file for file in os.listdir(f"uves_data/{folder}") if file.endswith("_bary.fits")]
            os.makedirs(f"spectra_UVES_bary/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:

                hdul = fits.open(f"uves_data/{folder}/{file}")
                data = hdul[1].data

                wavelength = data["WAVE"][0]
                flux = data["FLUX"][0]


                plt.plot(data["WAVE"][0], data["FLUX"][0], linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)")
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_UVES_bary/{folder}/{file}_spectra_UVES_bary.png")
                plt.close()
                print(f"UVES: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(data["FLUX"][0])


            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_uves.append((folder, common_wave, average_flux))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)")
            plt.savefig(f"combined_spectra_bary/UVES/average_spectrum_{folder}_bary.png")
            plt.close()
            print(f"UVES: Combined spectra {folder} saved.")



    if len(feros_array) == 0:
        print("No feros data found. Skipping")

    else:


        print(feros_array)
        for folder in feros_array:

            folder_data = [file for file in os.listdir(f"feros_data/{folder}") if file.endswith("_bary.fits")]
            os.makedirs(f"spectra_FEROS_bary/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:

                hdul = fits.open(f"feros_data/{folder}/{file}")
                data = hdul[1].data



                plt.plot(data["WAVE"][0], data["FLUX"][0], linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux (ADU)")
                plt.ylim(np.nanpercentile(data["FLUX"][0], [1, 99]))
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_FEROS_bary/{folder}/{file}_spectra_FEROS_bary.png")
                plt.close()
                print(f"FEROS: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(data["FLUX"][0])


            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_feros.append((folder, common_wave, average_flux))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux (ADU)")
            plt.ylim(np.nanpercentile(average_flux, [1, 99]))
            plt.savefig(f"combined_spectra_bary/FEROS/average_spectrum_{folder}_bary.png")
            plt.close()
            print(f"FEROS: Combined spectra {folder} saved.")





    flux_labels = {
        "UVES": "Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)",
        "FEROS": "Flux (ADU)"
    }

    wave_ranges = {
        "UVES": [("blue", 0, 5760), ("red", 5760, np.inf)],
        "FEROS": [("full", 0, np.inf)]
    }

    for instrument, grid_spectra in [
        ("UVES", grid_spectra_uves),
        ("FEROS", grid_spectra_feros)
    ]:

        for range_name, wave_lo, wave_hi in wave_ranges[instrument]:

            masks = [(wave >= wave_lo) & (wave < wave_hi) & (flux != 0) for _, wave, flux in grid_spectra]

            x_min = min(np.min(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            x_max = max(np.max(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_min = min(np.nanpercentile(flux[mask], 1) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_max = max(np.nanpercentile(flux[mask], 99) for (_, wave, flux), mask in zip(grid_spectra, masks))

            for start in range(0, len(grid_spectra), 8):

                group = grid_spectra[start:start + 8]
                nrows = (len(group) + 1) // 2

                fig, axes = plt.subplots(nrows, 2, figsize=(14, 3 * nrows), sharex=True, sharey=True, squeeze=False)
                axes = axes.flatten()

                for ax, (folder, wavelength, flux) in zip(axes, group):
                    ax.plot(wavelength, flux, linewidth=0.7)
                    ax.set_title(folder, fontsize=9)
                    ax.set_xlim(x_min, x_max)
                    ax.set_ylim(y_min, y_max)
                    ax.grid(True, alpha=0.2)

                for ax in axes[len(group):]:
                    ax.set_visible(False)

                fig.supxlabel("Wavelength (Å)")
                fig.supylabel(flux_labels[instrument])
                fig.tight_layout()

                fig.savefig(
                    f"combined_spectra_bary/{instrument}_spectral_grid_{range_name}_{start // 8 + 1}_bary.png",
                    dpi=300,
                    bbox_inches="tight"
                )
                plt.close(fig)
























def population():
    from astroquery.gaia import Gaia
    
    bg = Gaia.launch_job_async("""
    select top 100 parallax, phot_g_mean_mag, bp_rp
    from gaiadr3.gaia_source
    where parallax > 6.6667
      and parallax_over_error > 10
      and bp_rp is not null
      and phot_g_mean_mag is not null
    order by random_index
    """).get_results().to_pandas()

    ids = pd.read_csv('Gaia_Data_DR3.csv', dtype={'source_id': 'int64'})['source_id']
    cand = Gaia.launch_job(f"""
    select gs.source_id, gs.parallax, gs.phot_g_mean_mag, gs.bp_rp, bm.m1, bm.m1_lower, bm.m1_upper, bm.m2, bm.m2_lower, bm.m2_upper, nss.period, nss.semi_amplitude_primary, nss.significance, nss.goodness_of_fit
    from gaiadr3.gaia_source as gs
    join gaiadr3.binary_masses as bm
        on gs.source_id = bm.source_id
    join gaiadr3.nss_two_body_orbit as nss
        on gs.source_id = nss.source_id
    where gs.source_id in ({','.join(map(str, ids))})
    """).get_results().to_pandas()


    pop = Gaia.launch_job_async("""
    select bm.m1, bm.m2, bm.m2_lower, nss.period, nss.semi_amplitude_primary, nss.significance, nss.goodness_of_fit
    from gaiadr3.binary_masses as bm
    left join gaiadr3.nss_two_body_orbit as nss
        on bm.source_id = nss.source_id
    """).get_results().to_pandas()
    

    for df in (bg, cand):
        df.columns = df.columns.str.lower()
        df['MG'] = df['phot_g_mean_mag'] + 5 * np.log10(df['parallax']) - 10

    plt.figure(figsize=(6, 7))
    plt.scatter(bg['bp_rp'], bg['MG'], s=0.5, c='0.7', lw=0, rasterized=True, label='Gaia DR3 sample')
    uves = cand['source_id'].isin(pd.read_csv('UVES_data.csv')['GAIA_SOURCE_ID'])
    feros = cand['source_id'].isin(pd.read_csv('FEROS_data.csv')['GAIA_SOURCE_ID'])

    plt.scatter(cand['bp_rp'][uves & ~feros], cand['MG'][uves & ~feros], s=40, c='#1baf7a', edgecolor='k', zorder=3, label='UVES')
    plt.scatter(cand['bp_rp'][feros & ~uves], cand['MG'][feros & ~uves], s=40, c='#eb6834', edgecolor='k', zorder=3, label='FEROS')
    plt.scatter(cand['bp_rp'][uves & feros], cand['MG'][uves & feros], s=40, c='#2a78d6', edgecolor='k', zorder=3, label='Both')
    plt.gca().invert_yaxis()
    plt.xlabel('BP_RP')
    plt.ylabel('M_G')
    plt.legend()
    plt.savefig("candidates_on_CMD.png")
    

    plt.figure()
    plt.scatter(pop['m1'], pop['m2'], s=0.5, c='0.7', lw=0, rasterized=True, label='Gaia DR3 binary masses')
    plt.scatter(cand['m1'][uves & ~feros], cand['m2_lower'][uves & ~feros], s=40, c='#1baf7a', edgecolor='k', zorder=3, label='UVES')
    plt.scatter(cand['m1'][feros & ~uves], cand['m2_lower'][feros & ~uves], s=40, c='#eb6834', edgecolor='k', zorder=3, label='FEROS')
    plt.scatter(cand['m1'][uves & feros], cand['m2_lower'][uves & feros], s=40, c='#2a78d6', edgecolor='k', zorder=3, label='Both')

    plt.errorbar(cand['m1'][uves & ~feros], cand['m2_lower'][uves & ~feros], xerr=[cand['m1'][uves & ~feros] - cand['m1_lower'][uves & ~feros], cand['m1_upper'][uves & ~feros] - cand['m1'][uves & ~feros]], fmt='none', ecolor='k', capsize=3, zorder=2)
    plt.errorbar(cand['m1'][feros & ~uves], cand['m2_lower'][feros & ~uves], xerr=[cand['m1'][feros & ~uves] - cand['m1_lower'][feros & ~uves], cand['m1_upper'][feros & ~uves] - cand['m1'][feros & ~uves]], fmt='none', ecolor='k', capsize=3, zorder=2)
    plt.errorbar(cand['m1'][uves & feros], cand['m2_lower'][uves & feros], xerr=[cand['m1'][uves & feros] - cand['m1_lower'][uves & feros], cand['m1_upper'][uves & feros] - cand['m1'][uves & feros]], fmt='none', ecolor='k', capsize=3, zorder=2)
    plt.yscale('log')
    plt.xlabel('M1')
    plt.ylabel('M2')
    plt.legend()
    plt.savefig("candidates_m2_vs_m1.png")
    plt.close()


    plt.figure()
    plt.scatter(pop['period'], pop['semi_amplitude_primary'], s=0.5, c='0.7', lw=0, rasterized=True, label='Gaia DR3 binary masses')
    plt.scatter(cand['period'][uves & ~feros], cand['semi_amplitude_primary'][uves & ~feros], s=40, c='#1baf7a', edgecolor='k', zorder=3, label='UVES')
    plt.scatter(cand['period'][feros & ~uves], cand['semi_amplitude_primary'][feros & ~uves], s=40, c='#eb6834', edgecolor='k', zorder=3, label='FEROS')
    plt.scatter(cand['period'][uves & feros], cand['semi_amplitude_primary'][uves & feros], s=40, c='#2a78d6', edgecolor='k', zorder=3, label='Both')
    plt.yscale('log')
    plt.xlabel('Orbital Period')
    plt.ylabel('Velocity Semi-amplitude')
    plt.legend()
    plt.savefig("candidates_semi_amp_vs_period.png")
    plt.close()


    plt.figure()
    plt.scatter(pop['goodness_of_fit'], pop['significance'], s=0.5, c='0.7', lw=0, rasterized=True, label='Gaia DR3 binary masses')
    plt.scatter(cand['goodness_of_fit'][uves & ~feros], cand['significance'][uves & ~feros], s=40, c='#1baf7a', edgecolor='k', zorder=3, label='UVES')
    plt.scatter(cand['goodness_of_fit'][feros & ~uves], cand['significance'][feros & ~uves], s=40, c='#eb6834', edgecolor='k', zorder=3, label='FEROS')
    plt.scatter(cand['goodness_of_fit'][uves & feros], cand['significance'][uves & feros], s=40, c='#2a78d6', edgecolor='k', zorder=3, label='Both')
    plt.xlim(-5, 25)
    plt.yscale('log')
    plt.xlabel('Goodness of Fit')
    plt.ylabel('Significance')
    plt.legend()
    plt.savefig("candidates_signif_vs_goodness.png")
    plt.close()





def normalised_spectra():
    

    uves_array = os.listdir("uves_data")
    os.makedirs("spectra_UVES_normal", exist_ok=True)

    os.makedirs("combined_spectra_normal", exist_ok=True)
    os.makedirs("combined_spectra_normal/FEROS", exist_ok=True)
    os.makedirs("combined_spectra_normal/UVES", exist_ok=True)

    feros_array = os.listdir("feros_data")
    os.makedirs("spectra_FEROS_normal", exist_ok=True)

    grid_spectra_uves = []
    grid_spectra_feros = []

    if len(uves_array) == 0:
        print("No UVES data found. Skipping")

    else:


        print(uves_array)
        peak_rows = []
        for folder in uves_array:

            folder_data = [file for file in os.listdir(f"uves_data/{folder}") if file.endswith("_bary.fits")]
            os.makedirs(f"spectra_UVES_normal/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:

                hdul = fits.open(f"uves_data/{folder}/{file}")
                data = hdul[1].data

                wavelength = data["WAVE"][0]
                flux_base = data["FLUX"][0]

                flux = np.zeros_like(flux_base)

                for chip in [wavelength < 5760, wavelength >= 5760]:
                    use = chip & (flux_base != 0)
                    fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, niter=10, sigma_lower=1.5, sigma_upper=3)
                    continuum, _ = fitter(models.Chebyshev1D(10, domain=[wavelength[use].min(), wavelength[use].max()]), wavelength[use], flux_base[use])
                    flux[use] = flux_base[use] / continuum(wavelength[use])


                plt.plot(wavelength, flux, linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux normalised")
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_UVES_normal/{folder}/{file}_spectra_UVES_normal.png")
                plt.close()
                print(f"UVES: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(flux)


            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_uves.append((folder, common_wave, average_flux))

            search_flux = np.where(average_flux == 0, 1, average_flux)
            peaks, _ = find_peaks(1 - search_flux, prominence=0.1, width=3)
            peak_rows.append([folder] + list(np.round(common_wave[peaks], 2)))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.plot(common_wave[peaks], average_flux[peaks], "rv", markersize=4)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux normalised")
            plt.savefig(f"combined_spectra_normal/UVES/average_spectrum_{folder}_normal.png")
            plt.close()
            print(f"UVES: Combined spectra {folder} saved.")

        with open("UVES_peaks.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["GaiaID"] + [f"wavelength{i + 1}" for i in range(max(len(row) - 1 for row in peak_rows))])
            writer.writerows(peak_rows)



            

    if len(feros_array) == 0:
        print("No feros data found. Skipping")

    else:


        print(feros_array)
        for folder in feros_array:

            folder_data = [file for file in os.listdir(f"feros_data/{folder}") if file.endswith("_bary.fits")]
            os.makedirs(f"spectra_FEROS_normal/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:

                hdul = fits.open(f"feros_data/{folder}/{file}")
                data = hdul[1].data
                wavelength = data["WAVE"][0]
                flux_base = data["FLUX"][0]

                flux = np.zeros_like(flux_base)

                for chip in [wavelength > 3800]:
                    use = chip & (flux_base != 0)
                    continuum = np.polyval(np.polyfit(wavelength[use], flux_base[use], 5), wavelength[use])
                    above = flux_base[use] > continuum
                    continuum = np.polyval(np.polyfit(wavelength[use][above], flux_base[use][above], 5), wavelength[use])
                    flux[use] = flux_base[use] / continuum



                plt.plot(data["WAVE"][0], flux, linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux normalised")
                plt.ylim(np.nanpercentile(flux, [1, 99]))
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_FEROS_normal/{folder}/{file}_spectra_FEROS_normal.png")
                plt.close()
                print(f"FEROS: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(flux)


            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_feros.append((folder, common_wave, average_flux))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux normalised")
            plt.ylim(np.nanpercentile(average_flux, [1, 99]))
            plt.savefig(f"combined_spectra_normal/FEROS/average_spectrum_{folder}_normal.png")
            plt.close()
            print(f"FEROS: Combined spectra {folder} saved.")





    flux_labels = {
        "UVES": "Flux normalised",
        "FEROS": "Flux normalised"
    }

    wave_ranges = {
        "UVES": [("blue", 0, 5760), ("red", 5760, np.inf)],
        "FEROS": [("full", 0, np.inf)]
    }

    for instrument, grid_spectra in [
        ("UVES", grid_spectra_uves),
        ("FEROS", grid_spectra_feros)
    ]:

        for range_name, wave_lo, wave_hi in wave_ranges[instrument]:

            masks = [(wave >= wave_lo) & (wave < wave_hi) & (flux != 0) for _, wave, flux in grid_spectra]

            x_min = min(np.min(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            x_max = max(np.max(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_min = min(np.nanpercentile(flux[mask], 1) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_max = max(np.nanpercentile(flux[mask], 99) for (_, wave, flux), mask in zip(grid_spectra, masks))

            for start in range(0, len(grid_spectra), 8):

                group = grid_spectra[start:start + 8]
                nrows = (len(group) + 1) // 2

                fig, axes = plt.subplots(nrows, 2, figsize=(14, 3 * nrows), sharex=True, sharey=True, squeeze=False)
                axes = axes.flatten()

                for ax, (folder, wavelength, flux) in zip(axes, group):
                    ax.plot(wavelength, flux, linewidth=0.7)
                    ax.set_title(folder, fontsize=9)
                    ax.set_xlim(x_min, x_max)
                    ax.set_ylim(0, 1.4)
                    ax.grid(True, alpha=0.2)

                for ax in axes[len(group):]:
                    ax.set_visible(False)

                fig.supxlabel("Wavelength (Å)")
                fig.supylabel(flux_labels[instrument])
                fig.tight_layout()

                fig.savefig(
                    f"combined_spectra_normal/{instrument}_spectral_grid_{range_name}_{start // 8 + 1}_normal.png",
                    dpi=300,
                    bbox_inches="tight"
                )
                plt.close(fig)







def inspector():


    target = "4973506334166805376"

    folder_data = os.listdir(f"uves_data/{target}")

    for file in folder_data:

        if file.endswith("_bary.fits"):

            hdul = fits.open(f"uves_data/{target}/{file}")
            data = hdul[1].data

            wavelength = data["WAVE"][0]
            flux_base = data["FLUX"][0]

            flux = np.zeros_like(flux_base)

            for chip in [wavelength < 5760, wavelength >= 5760]:
                use = chip & (flux_base != 0)
                fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, niter=10, sigma_lower=1.5, sigma_upper=3)
                continuum, _ = fitter(models.Chebyshev1D(10, domain=[wavelength[use].min(), wavelength[use].max()]), wavelength[use], flux_base[use])
                flux[use] = flux_base[use] / continuum(wavelength[use])

            plt.plot(wavelength, flux, linewidth=1, label=file)

        else:
            continue

    plt.xlim(5420, 5430)
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux normalised")
    plt.savefig("4973506334166805376.png")
    plt.show()


def spectral_lines():

    target = "4239611700214875904"

    folder_data = os.listdir(f"uves_data/{target}")

    peak_rows = []
    spectra = []
    wavelengths = []

    for file in folder_data:
        if file.endswith("_bary.fits"):
                
            folder = target
            hdul = fits.open(f"uves_data/{folder}/{file}")
            data = hdul[1].data

            wavelength = data["WAVE"][0]
            flux_base = data["FLUX"][0]

            flux = np.zeros_like(flux_base)

            for chip in [wavelength < 5760, wavelength >= 5760]:
                use = chip & (flux_base != 0)
                fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, niter=10, sigma_lower=1.5, sigma_upper=3)
                continuum, _ = fitter(models.Chebyshev1D(10, domain=[wavelength[use].min(), wavelength[use].max()]), wavelength[use], flux_base[use])
                flux[use] = flux_base[use] / continuum(wavelength[use])

            wavelengths.append(data["WAVE"][0])
            spectra.append(flux)


    common_wave = wavelengths[0]
    spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
    average_flux = np.mean(spectra, axis=0)

    search_flux = np.where(average_flux == 0, 1, average_flux)
    peaks, _ = find_peaks(1 - search_flux, prominence=0.1, width=3)
    peak_rows.append([folder] + list(np.round(common_wave[peaks], 2)))

    plt.plot(common_wave, average_flux, linewidth=1)





    lines = {
        4861.34: r'F (H$\beta$)',
        4957.61: 'c (Fe)',
        5167.33: r'b$_4$ (Mg)',
        5168.91: r'b$_3$ (Fe)',
        5172.70: r'b$_2$ (Mg)',
        5183.62: r'b$_1$ (Mg)',
        5270.39: r'E$_2$ (Fe)',
        5460.73: 'e (Hg)'
    }
    
    for wavelength, label in lines.items():
        plt.axvline(x=wavelength, color='red', alpha=0.25, linestyle='--')
        plt.text(wavelength, 0.05, label, rotation=90, color='red',
                 transform=plt.gca().get_xaxis_transform(), fontsize=8,
                 va='bottom', ha='right')







    plt.xlim(4800, 5720)
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux normalised")
    plt.savefig("spectral_lines_pre.png")
    plt.show()







    peak_rows = []
    spectra = []
    wavelengths = []

    for file in folder_data:
        folder = target
        hdul = fits.open(f"uves_data/{folder}/{file}")
        data = hdul[1].data

        wavelength = data["WAVE"][0]
        flux_base = data["FLUX"][0]

        flux = np.zeros_like(flux_base)

        for chip in [wavelength < 5760, wavelength >= 5760]:
            use = chip & (flux_base != 0)
            fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, niter=10, sigma_lower=1.5, sigma_upper=3)
            continuum, _ = fitter(models.Chebyshev1D(10, domain=[wavelength[use].min(), wavelength[use].max()]), wavelength[use], flux_base[use])
            flux[use] = flux_base[use] / continuum(wavelength[use])

        wavelengths.append(data["WAVE"][0])
        spectra.append(flux)


    common_wave = wavelengths[0]
    spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
    average_flux = np.mean(spectra, axis=0)

    search_flux = np.where(average_flux == 0, 1, average_flux)
    peaks, _ = find_peaks(1 - search_flux, prominence=0.1, width=3)
    peak_rows.append([folder] + list(np.round(common_wave[peaks], 2)))

    plt.plot(common_wave, average_flux, linewidth=1)





    lines = {
        5875.618: r'D$_3$ (He)',
        5889.95: r'D$_2$ (Na)',
        5895.92: r'D$_1$ (Na)',
        6276.61: r'a (O$_2$)',
        6562.81: r'C (H$\alpha$)'
    }

    for wavelength, label in lines.items():
        plt.axvline(x=wavelength, color='red', alpha=0.25, linestyle='--')
        plt.text(wavelength, 0.05, label, rotation=90, color='red',
                 transform=plt.gca().get_xaxis_transform(), fontsize=8,
                 va='bottom', ha='right')








    plt.xlim(5820, 6780)
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux normalised")
    plt.savefig("spectral_lines_post.png")
    plt.show()


def spectrum_splitter():

    uves_array = os.listdir("uves_data")


    for folder in uves_array:

        folder_data = os.listdir(f"uves_data/{folder}")

        for file in folder_data:

            if file.endswith("_bary.fits"):

                chip_lower = 5750
                chip_upper = 5850

                hdul = fits.open(f"uves_data/{folder}/{file}")
                data = Table.read(f"uves_data/{folder}/{file}", hdu=1)
        

                wavelength = data["WAVE"]

                uves_blue = data[wavelength < chip_lower]
                uves_red = data[wavelength > chip_upper]

                uves_blue.write(f"uves_data/{folder}/{file[:-5]}_blue.fits", overwrite=True)
                uves_red.write(f"uves_data/{folder}/{file[:-5]}_red.fits", overwrite=True)





def error_pdf():

    uves_array = os.listdir("uves_data")

    for target in uves_array:

        folder_data = os.listdir(f"uves_data/{target}")

        flat_band_lower_blue = 5400
        flat_band_higher_blue = 5450

        flat_band_lower_red = 6050
        flat_band_higher_red = 6100


        res_red_array = []
        res_blue_array = []

        pix_blue_array = []
        pix_red_array = []
        
        for file in folder_data:
            if file.endswith("_bary.fits"):
                    
                folder = target
                hdul = fits.open(f"uves_data/{folder}/{file}")
                data = hdul[1].data


                wavelength = data["WAVE"][0]
                flux = data["FLUX"][0]

                err = data["ERR"][0]

                flat_band_blue = (wavelength > flat_band_lower_blue) & (wavelength < flat_band_higher_blue)
                flat_band_red = (wavelength > flat_band_lower_red) & (wavelength < flat_band_higher_red)

                fit_blue = np.polyval(np.polyfit(wavelength[flat_band_blue], flux[flat_band_blue], 1), wavelength[flat_band_blue])
                fit_red = np.polyval(np.polyfit(wavelength[flat_band_red], flux[flat_band_red], 1), wavelength[flat_band_red])
                residual_blue = (flux[flat_band_blue] - fit_blue) / err[flat_band_blue]
                residual_red = (flux[flat_band_red] - fit_red) / err[flat_band_red]
                #print(f"{target}, {file}, pixels: {np.sum(flat_band_blue)}, zero err: {np.sum(err[flat_band_blue] == 0)}, residual blue: {np.std(residual_blue)}")
                #print(f"Pixel difference blue: {np.std(np.diff(flux[flat_band_blue])) / np.sqrt(2) / np.median(err[flat_band_blue])}")
                #print(f"Pixel difference red: {np.std(np.diff(flux[flat_band_red])) / np.sqrt(2) / np.median(err[flat_band_red])}")


                plt.hist(np.diff(flux[flat_band_blue]) / (np.sqrt(2) * err[flat_band_blue][:-1]), bins = 30, density = True, alpha = 0.5, label = "Blue chip")
                plt.hist(np.diff(flux[flat_band_red]) / (np.sqrt(2) * err[flat_band_red][:-1]), bins = 30, density = True, alpha = 0.5, label = "Red chip")
                
                x = np.linspace(-5, 5, 500)

                plt.plot(x, np.exp(-x**2 / 2) / np.sqrt(2 * np.pi), label = "N(0, 1)")
                sigma_blue = np.std(np.diff(flux[flat_band_blue]) / (np.sqrt(2) * err[flat_band_blue][:-1]))
                sigma_red = np.std(np.diff(flux[flat_band_red]) / (np.sqrt(2) * err[flat_band_red][:-1]))
                plt.plot(x, np.exp(-x**2 / (2 * sigma_blue**2)) / (sigma_blue * np.sqrt(2 * np.pi)), linestyle = "--", label = f"N(0, {sigma_blue:.2f}) blue")
                plt.plot(x, np.exp(-x**2 / (2 * sigma_red**2)) / (sigma_red * np.sqrt(2 * np.pi)), linestyle = ":", label = f"N(0, {sigma_red:.2f}) red")
                plt.xlabel(r"$\Delta F / (\sqrt{2}\,\sigma)$")
                plt.ylabel("Probability density")
                plt.legend()
                #plt.show()

                res_blue_array.append(np.std(residual_blue))
                res_red_array.append(np.std(residual_red))

                pix_blue_array.append(sigma_blue)
                pix_red_array.append(sigma_red)
                
                #print(f"Residual blue: {np.std(residual_blue)}, {folder}, {file}")
                #print(f"Residual red: {np.std(residual_red)}, {folder}, {file}")
       
        print(f"Target: {target}: Residual blue mean: {np.mean(res_blue_array)}, Residual blue std: {np.std(res_blue_array)}")
        print(f"Target: {target}: Residual red mean: {np.mean(res_red_array)}, Residual red std: {np.std(res_red_array)}")


        print(f"Target: {target}: Pixel difference blue mean: {np.mean(pix_blue_array)}, Pixel difference blue std: {np.std(pix_blue_array)}")
        print(f"Target: {target}: Pixel difference red mean: {np.mean(pix_red_array)}, Pixel difference red std: {np.std(pix_red_array)}")

#downloading_data()
#population()
#barycentric_corrections()
#plot_spectra_old()
#plot_spectra_bary()
#normalised_spectra()
inspector()
#spectral_lines()
#spectrum_splitter()
#error_pdf()

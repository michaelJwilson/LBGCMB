;; name: synthesise_mags.pro
;;
;; aim: to synthesise the mags of the dwarf stars for injection into
;; my simulations.
;;
;; created: Thursday 8th May 2014
;;
pro synthesise_mags

filesave = 'dwarfStar_mags.sav'

;; read in the required filters
filtersFile = 'filtersFile.lis'
readcol, filtersFile, filterNames, filterTransmissionName, format = 'A, A', skipline = 1
numberFilters = n_elements(filterNames)

filtersDir = '/disk2/raab/photo_z/filters/'

;; read in the stellar spectra
starDictionaryFile = '/disk2/raab/photo_z/lePhare/lephare_dev/sed/STAR/DWARFSTARS/DWARFSTARS_MOD.list'

;; read in the star list
readcol, starDictionaryFile, starFilename,ascii_type, format = "A, A", count = numberStarTemplates

directory = '/disk1/raab/disk3/lePhare/lephare_dev/sed/STAR/'

;; extract the star type from the name
starType = strmid(starFilename, 11, 2)
starIDs  = indgen(n_elements(starFilename))+1

starMags = fltarr(numberStarTemplates, numberFilters)

;; loop through each star
for si = 0, numberStarTemplates -1 do begin
   print, "Extracting the mags for stellar type ", starType[si]
   
;; read in the star spectrum and synthesise:
   print, directory + starFilename[si]
   readcol, directory + starFilename[si], wavelength, flux_lambda, $
            format = 'F, F'
   
   cgplot, wavelength, flux_lambda, xrange = [0, 20000]
   
;; loop through the filters necessary, integrating through
   for fi = 0, numberFilters -1 do begin
      print, "Synthesising mag for filter: ", filterNames[fi]
      
      starMags[si, fi] = observethroughfilter(wavelength, flux_lambda, filterTransmissionName[fi], /returnflux)
      
      if starMags[si,fi] lt 0.0 then starMags[si,fi] = 0.0
      print,  starMags[si, fi]
      
      
      
   endfor
   
  ;; cafe_pause
   
   
   
endfor

save, starMags, starType, filterNames, filename = filesave
print, "Results saved to ", filesave


end

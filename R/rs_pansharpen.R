rs_pansharpen <- function(bands,
                          panchro,
                          filter = matrix(c(1,1,1,1,2,1,1,1,1),
                                          nrow = 3),
                          scale = 2){

  # Resample ----------------------------------
  resampled <- purrr::map(bands, terra::resample, panchro)

  # Matching ----------------------------------
  # mean and sd of panchromatic
  mean_pan <- terra::global(panchro, mean)[,1]
  sd_pan <- terra::global(panchro, sd)[,1]

  # Map matching function
  matched_pan <- purrr::map(resampled, function(band){

    # mean and sd of band i
    mean_band <- terra::global(band, mean)[,1]
    sd_band <- terra::global(band, sd)[,1]

    # Eq2
    ai <- sd_band / sd_pan
    # Eq3
    bi <- mean_band - (ai * mean_pan)
    # Eq1
    pan_xi <- ai * panchro + bi

    # Result of matching
    return(pan_xi)
  })

  # Filter
  for (i in 2:scale){
    if (i == 2){
      filtered_pan <- purrr::map(matched_pan,
                          function(band){
                            lpf <- terra::focal(band,
                                         filter = filter,
                                         fun = mean)
                            return(lpf)})
    } else {
      filtered_pan <- purrr::map(filtered_pan,
                          function(band){
                            lpf <- terra::focal(band,
                                         filter = filter,
                                         fun = mean)
                            return(lpf)})
    }}

  # Detail
  detail <- purrr::map2(matched_pan,
                 filtered_pan,
                 function(matchBand, filterBand){
                   detail <- matchBand - filterBand
                   return(detail)
                 })
  # Fusion
  fusion <- purrr::map2(resampled,
                        detail,
                        function(x, y){return(x + y)})

  return(fusion)
}

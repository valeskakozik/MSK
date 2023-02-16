EDF_extr_EEG = function(d){
  
  options(scipen = 999) #disable scientific notation
  
  cat('\f')
  
if ( !require(edf) ) {cat('\n'); install.packages("edf")}
if ( !require(edfReader) ) {cat('\n'); install.packages("edfReader")}
if ( !require(purrr) ) {cat('\n'); install.packages("purrr")}
  
  cat('\n', paste(rep('= ', 30), collapse = ""),
      '\n\n\t\t\tEDF_extr_EEG(d)\n',
      '\n\t\tValeska Kozik\t\t21-04-2021\n\n',
      paste(rep('= ', 30), collapse = ""),
      '\n\nThis script extracts EEG signal and metadata from .edf files.',
      '\nDo not close this R session until the script is done.\n')

#...Setup directory
setup = list(
  directory = dirname(rstudioapi::getActiveDocumentContext()$path)
  )
#setup = list(directory = d)
setup$source = paste0(setup$directory, .Platform$file.sep, 'EDF_FILES')
  
#...Create destination directories
setup$dest = list(info   = paste0(setup$directory, .Platform$file.sep, '00_info'),
                  annot  = paste0(setup$directory, .Platform$file.sep, '01_annotations'),
                  signal = paste0(setup$directory, .Platform$file.sep, '02_signal'),
                  temp   = paste0(setup$directory, .Platform$file.sep, '_out_temp'))

if (length(setdiff(setup$dest, list.dirs(setup$directory))) != 0 ) {
  cat('\n\nCreating directories ...')
  invisible(lapply(setdiff(setup$dest, list.dirs(setup$directory)), dir.create)); cat(' Done.\n')
  }

#...List files that match pattern *.edf
setup$edf.files = dir(setup$source,'*.EDF')
cat(paste0('\nPattern: *.edf\nFiles matching pattern: ', length(setup$edf.files), '\n'))

setup$electrodes = matrix(nrow = length(setup$edf.files), ncol = 2)
electrodes = c('Cz', 'F3', 'F4', 'P3', 'P4', 'T3', 'T4', 'O1', 'O2')

for (i in 1:length(setup$edf.files)){
  
  #... Read in EDF-file
  cat(paste0('\nLoading file [', i, '/', length(setup$edf.files), ']: ',  setup$edf.files[i], ' ...'))            
  file = edf::read.edf(paste(setup$source, setup$edf.files[i], sep =.Platform$file.sep))
  header = edfReader::readEdfHeader(paste(setup$source, setup$edf.files[i], sep =.Platform$file.sep))
  
  #...Change electrode name from "EEG1" to "Cz"
  names(file$header.signal)[grep('EEG1', names(file$header.signal))] = 
  names(file$signal)[grep('EEG1', names(file$signal))]               = 'Cz'
  
  #...carry over only relevant electrodes
  file$signal = file$signal[match(electrodes, names(file$signal))]
  file$header.signal = file$header.signal[match(electrodes, names(file$header.signal))]
  
  #...Fetch metadata and annotations per file
  meta = list( 
    rec.ID        = file$header.global$patient.id,
    pat.ID        = regmatches(file$header.global$patient.id, regexpr('(?=MSK).+', file$header.global$patient.id, perl = T)),
    rec.date      = format(file$header.global$timestamp.start, '%d %b %Y'),
    rec.time_secs = as.numeric(difftime(file$header.global$timestamp.stop, file$header.global$timestamp.start, units='secs')),
    freq          = file$header.signal[[1]]$samplingrate,
    dim           = file$header.signal[[1]]$physical.dimension,
    n.annot       = length(file$events[file$events$annotation != "", "onset"])
    )
  
  setup$electrodes[i, 1:2] = c(meta$pat.ID, paste(names(file$signal), collapse = ' - '))
  
  prefiltering  = data.frame( t( sapply(file$header.signal, `[[`, 8) ) )
  names(prefiltering) = paste0('prefilter.', names(prefiltering))
  
  if ( i == 1 ){
    meta.table = cbind(data.frame(meta), prefiltering)
  }else{
    meta.table = rbind(meta.table, cbind(data.frame(meta), prefiltering))
  }
  
  #...Fetch signal for all electrodes and write to ASCII file
  cat('\tcreating output files ...')
  signal = data.frame(sapply(file$signal, `[[`, 1))

  write.table(
    signal,
    paste0(setup$dest$temp, .Platform$file.sep, meta$pat.ID, '_signal.asc'),
    dec = '.',
    sep = '\t',
    row.names = F,
    col.names = F
  )
  
  #... Write annotations to files
  write.csv(
    data.frame(onset_sec = file$events[file$events$annotation != "", "onset"],
               label = file$events$annotation[file$events$annotation != ""]),
    paste0(setup$dest$annot, .Platform$file.sep, meta$pat.ID, '_annotations.txt'),
    row.names = F
    ); cat('\tDone.')

#... Write metadata to file
  if ( i == length(setup$edf.files) ){
    cat('\nCreating metadata file ... ')
    write.csv(
      meta.table,
      paste0(setup$dest$info, .Platform$file.sep, 'EEG_metadata_edf_', format(Sys.Date(), '%d_%m_%Y'), '.csv'),
      row.names = F
    )
    
    if ( length(unique(meta.table$electrodes)) == 1 ){
      write.table(
        t(matrix(names(file$signal))),
        paste(setup$dest$info, 'electrodes.txt', sep = .Platform$file.sep),
        sep = ',',
        row.names = F,
        col.names = F
      )
      }else{
       #write unique to table 
        elc = matrix(nrow = length(unique(setup$electrodes[,2])), ncol = 2)
        for ( j in 1 : length(unique(meta.table$electrodes)) ){
          elc[j, 1] = paste(setup$electrodes[, 1][setup$electrodes == unique(setup$electrodes)[j, 1]], collapse = '; ')
          elc[j, 2] = unique(setup$electrodes[, 2])[j]
        }
        write.table(
          elc,
          paste(setup$dest$info, 'electrodes.txt', sep = .Platform$file.sep),
          sep = ',',
          row.names = F,
          col.names = F
        )
      }
    
    cat('Done.\n\nCleaning up ...')
    
    if ( length(setup$edf.files) == length(dir(setup$dest$temp, '*.asc')) ){
      
      invisible(file.copy(
        paste(setup$dest$temp, dir(setup$dest$temp), sep = .Platform$file.sep),
        paste(setup$dest$signal, dir(setup$dest$temp), sep = .Platform$file.sep)
      )); cat('\tDone.\n')
      
      cat(paste0('\n[', length(dir(setup$dest$temp)), '/', length(setup$edf.files), '] files converted.\n',
      '\nOutput locations:', '\nMetadata:\t', setup$dest$info, '\nAnnotations:\t', setup$dest$annot, '\nASCII files:\t', setup$dest$signal))
      
      cat('\n\nDone! You may shut down R.\n')
      
      unlink(setup$dest$temp, recursive = T)
      
      o.warn = getOption('warn')
      options(warn = -1)
    }else{
      cat('\nError: Export unsuccessful!')
    }
    }
}
options(scipen = 0) #turn scientific notation back on

on.exit(options(warn = o.warn))
}

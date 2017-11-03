#' @export convert_raw_fid_to_complex
convert_raw_fid_to_complex <- function(rawValue, endian="little") {
  
  lengthraw <- length(rawValue)
  fid <- readBin(rawValue, 'double', n = length(rawValue)/4, size = 4, signed = TRUE)
  odd <- seq(1, lengthraw/4, by=2)
  even <- seq(2, lengthraw/4, by=2)
  complex(real=fid[odd], imaginary=fid[even])
}

parse_siemens_mrs_dicom_header <- function(rawString, sq.txt="", endian="little",
                             verbose=FALSE) {
  ##
  ## "The default DICOM Transfer Syntax, which shall be supported by
  ## all AEs, uses Little Endian encoding and is specified in Annex
  ## A.1." (PS 3.5-2004, page 38)
  ##
  ## PS 3.5-2004, Sect 7.1.2: Data Element Structure with Explicit VR
  ## Explicit VRs store VR as text chars in 2 bytes.
  ## VRs of OB, OW, SQ, UN, UT have VR chars, then 0x0000, then 32 bit VL:
  ##
  ## +-----------------------------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 | 10 | 11 |
  ## +----+----+----+----+----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<VR----->|<0x0000->|<Length----------->|<Value->
  ##
  ## Other Explicit VRs have VR chars, then 16 bit VL:
  ##
  ## +---------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
  ## +----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<VR----->|<Length->|<Value->
  ##
  ## Implicit VRs have no VR field, then 32 bit VL:
  ##
  ## +---------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
  ## +----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<Length----------->|<Value->
  ##
  rawToHex <- function(bytes) {
    toupper(paste(rev(bytes), collapse=""))
  }
  is.item <- function (group, element) {
    group == "FFFE" && element %in% c("E000","E00D","E0DD")
  }
  ## Data files that are necessary to proceed
  #data(dicom.dic, package="oro.dicom", envir=environment())
  #data(dicom.VR, package="oro.dicom", envir=environment())
  ##
  strseek <- dseek <- 0
  dicomHeader <- NULL
  pixelData <- spectroscopyData <- FALSE
  while (! (pixelData || spectroscopyData) && strseek < length(rawString)) {
    ## rm(group, element, dictionaryIndex, dic, rawValue, VR, vr, value, length)
    group <- rawToHex(rawString[strseek + 1:2])
    element <- rawToHex(rawString[strseek + 3:4])
    if (! any(dictionaryIndex <- group == oro.dicom::dicom.dic$group & element == oro.dicom::dicom.dic$element)) {
      ## Private tag = Unknown
      dic <- list(group=group, element=element, code="UN", offset=1, name="Unknown")
    } else {
      dic <- oro.dicom::dicom.dic[dictionaryIndex, ]
    }
    if (verbose) {
      cat("#", group, element, dic$name, dic$code, sep="\t")
    }
    b56 <- .rawToCharWithEmbeddedNuls(rawString[strseek + 5:6])
    b78 <- readBin(rawString[strseek + 7:8], "integer", size=2, endian=endian)
    b47 <- readBin(rawString[strseek + 5:8], "integer", size=4, endian=endian, signed = TRUE)
    strseek <- strseek + 8
    if (b56 %in% c("OB","OW","SQ","UN","UT") && ! is.item(group, element)) {
      ## Explicit VR
      lengthraw <- readBin(rawString[strseek + 1:4], "integer", size=4, endian=endian, signed = TRUE)
      strseek <- strseek + 4
      if (b56 != "SQ") {
        rawValue <- rawString[strseek + 1:lengthraw]
      }
      vr <- b56
    } else {
      if (b56 %in% c("AE","AS","AT","CS","DA","DS","DT","FL","FD","IS","LO",
                     "LT","OF","PN","SH","SL","SS","ST","TM","UI","UL","US")
          && ! is.item(group, element)) {
        ## Explicit VR
        lengthraw <- b78
        if (dic$code != "SQ") {
          rawValue <- rawString[strseek + 1:lengthraw]
        }
        vr <- b56
      } else {
        ## Implicit VR
        lengthraw <- b47
        vr <- NULL
        if (is.item(group, element)) {
          lengthraw <- 0
          rawValue <- raw(0)
        } else {
          if (dic$code != "SQ") {
            rawValue <- rawString[strseek + 1:lengthraw]
          } else {
            vr <- "SQ"
          }
        }
      }
    }
    if (! is.null(vr)) {
      VR <- oro.dicom::dicom.VR[oro.dicom::dicom.VR$code == vr, ]
    } else if (dic$code != "UN") {
      VR <- oro.dicom::dicom.VR[oro.dicom::dicom.VR$code == dic$code, ]
    } else {
      VR <- oro.dicom::dicom.VR[oro.dicom::dicom.VR$code == "UN", ]
    }
    if (verbose) {
      cat("", VR$code, lengthraw, sep="\t")
    }
    if (sq.txt == "" && group == "7FE0" && element == "0010") {
      ## PixelData
      value <- "PixelData"
      pixelData <- TRUE
      dseek <- strseek
    } else {
      if (sq.txt == "" && group == "7FE1" && element == "1010") {
        ## SpectroscopyData
        value <- "SpectroscopyData"
        spectroscopyData <- TRUE
        dseek <- strseek
        # dseek <- strseek + 4 # HACK: not sure why I need to skip an extra four bytes
      } else {
        if (VR$code %in% c("UL","US")) { # (VR$code == "UL" || VR$code == "US") {
          value <- readBin(rawValue, "integer", n=lengthraw/VR$bytes,
                           size=VR$bytes, signed=FALSE, endian=endian)
        } else if (VR$code %in% c("SL","SS")) { # (VR$code == "SL" || VR$code == "SS") {
          value <- readBin(rawValue, "integer", n=lengthraw/VR$bytes,
                           size=VR$bytes, signed=TRUE, endian=endian)
        } else if (VR$code %in% c("FD","FL")) { # (VR$code == "FD" || VR$code == "FL") {
          value <- readBin(rawValue, "numeric", n=lengthraw/VR$bytes,
                           size=VR$bytes, signed=TRUE, endian=endian)
        } else if (VR$code %in% c("OB","OW")) { # (VR$code == "OB" || VR$code == "OW") {
          value <- .rawToCharWithEmbeddedNuls(rawValue)
        } else if (VR$code == "SQ") {
          value <- "Sequence"
        } else {
          if (lengthraw > 0) {
            tmpString <- .rawToCharWithEmbeddedNuls(rawValue)
            tmpString <- sub(" +$", "", tmpString)     # remove white space at end
            tmpString <- gsub("[\\/]", " ", tmpString) # remove "/"s
            tmpString <- gsub("[\\^]", " ", tmpString) # remove "^"s
          } else {
            tmpString <- ""
          }
          value <- tmpString
        }
      }
    }
    if (verbose) {
      cat("", value, sq.txt, sep="\t", fill=TRUE)
    }
    dicomHeaderRow <- c(group, element, dic$name, dic$code, lengthraw, value, sq.txt)
    dicomHeader <- rbind(dicomHeader, dicomHeaderRow)
    if (group == "FFFE" && element == "E0DD") { # SequenceDelimitationItem
      dseek <- strseek
      break
    }
    if (VR$code == "SQ") {
      groupElement <- paste("(", group, ",", element, ")", sep="")
      if (lengthraw > 0) {
        ## Pass length of bytes provided explicitly by the sequence tag
        dcm <- parse_siemens_mrs_dicom_header(rawString[strseek + 1:lengthraw],
                                paste(sq.txt, groupElement),
                                verbose=verbose)
      } else {
        ## Pass remaining bytes and look for SequenceDelimitationItem tag
        dcm <- parse_siemens_mrs_dicom_header(rawString[(strseek + 1):length(rawString)],
                                paste(sq.txt, groupElement),
                                verbose=verbose)
        lengthraw <- dcm$data.seek
      }
      dicomHeader <- rbind(dicomHeader, dcm$header)
    }
    strseek <- strseek + ifelse(lengthraw >= 0, lengthraw, 0)
  }
  list(header = dicomHeader, pixel.data = pixelData, data.seek = dseek,
       spectroscopy.data = spectroscopyData, rawfid = rawValue)
}

.rawToCharWithEmbeddedNuls <- function(str.raw, to="UTF-8") {
  iconv(rawToChar(str.raw[str.raw != as.raw(0)]), to=to)
}

#' @export
read_siemens_mrs_dicom <- function(fname, endian="little", boffset=NULL, debug=FALSE) {
  
  fsize <- file.info(fname)$size
  fd <- file(fname, 'rb')
  # fraw <- readBin(fname, "raw", n=as.integer(fsize), endian=endian)
  fraw <- readBin(fd, "raw", n=as.integer(fsize), endian=endian)
  if (is.null(boffset) && any(as.integer(fraw[1:128]) != 0)) {
    stop("Non-zero bytes are present in the first 128, please use\nboffset to skip the necessary number of bytes.")
  }
  skip128 <- fraw[1:128]
  skipFirst128 <- ifelse(any(as.logical(skip128)), FALSE, TRUE)
  if (debug) {
    cat("# First 128 bytes of DICOM header =", fill=TRUE)
    print(skip128)
  }
  DICM <- .rawToCharWithEmbeddedNuls(fraw[129:132]) == "DICM"
  if (debug) {
    cat("# DICM =", DICM, fill=TRUE)
  }
  #if (DICM) {
  #  if (.rawToCharWithEmbeddedNuls(fraw[129:132]) != "DICM") {
  #    stop("DICM != DICM")
  #  }
  #}
  dicomHeader <- sequence <- NULL
  seq.txt <- ""
  if (is.null(boffset)) {
    bstart <- 1 + ifelse(skipFirst128, 128, 0) + ifelse(DICM, 4, 0) # number of bytes to start
  } else {
    bstart <- boffset + 1
  }  
  dcm <- parse_siemens_mrs_dicom_header(fraw[bstart:fsize], seq.txt, endian=endian, verbose=debug)
  complex_fid <- convert_raw_fid_to_complex(dcm$rawfid)
  close(fd)
  return(list(fid = complex_fid, dcm = dcm))

}

#' @export
replace_fid_in_siemens_mrs_dicom <- function(filename, complex_fid) {
  
  if (! file.exists(filename))
    stop(sprintf('Dicom file %s does not exist', filename), call. = FALSE)
  
  mrs_dcm <- read_siemens_mrs_dicom(filename)
  old_real_part <- Re(mrs_dcm$fid)
  old_imag_part <- Im(mrs_dcm$fid)
  new_real_part <- Re(complex_fid)
  new_imag_part <- Im(complex_fid)
  
  # Check if the length of the new complex_fid is same as the length of the existing complex_fid  
  if ( (length(old_real_part) != length(new_real_part)) ||  (length(old_imag_part) != length(new_imag_part)) ) 
    stop(sprintf('The length of Spectroscopy fid signal does not match the signal length in file\n'), call. = FALSE)

  fd <- file(filename, 'r+b')
  # Skip to the dicom fid signal
  seek(fd, where = mrs_dcm$dcm$data.seek + 132, rw = "write")
  odd <- seq(1, length(complex_fid)*2, by=2)
  even <- seq(2, length(complex_fid)*2, by=2)
  y = numeric(length(complex_fid)*2)
  y[odd] <- new_real_part
  y[even] <- new_imag_part
  writeBin(y, fd, size = 4)
  close(fd)
}


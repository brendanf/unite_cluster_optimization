#!/bin/awk -f
# The first file to look through is the protax reference file.
# Its headerss are sh_accno, so use _ as field separator.
# In some cases accno also contains an underscore.
BEGIN{FS="_"}

# The second file is the Unite + INSD fasta release to look through for matching
# sequences.
# Its headers are accno|taxonomy|sh, so use | as field separator starting
# after the end of the (first) file.
ENDFILE{FS="[|]"}

# Process fasta header lines
/^>/ {
  if (FNR==NR) {
    # the first file

    # The accno is the second field.
    # When we search in the second file, it will also include the ">" at
    # the beginning of the header.
    accno=">" $2
    # If the accno contains an underscore, then it will be split into two fields.
    # Recombine them.
    if ($3) accno=accno "_" $3

    # protax has some duplicate sequences (same accno, but different SH version
    # numbers)Keep track of these
    if (tbd[accno]) {
      ndups++
    } else {
      # Add the accno to the list of sequences which are needed.
      tbd[accno]=1
      # Add the full header (to match the Protax file)
      header[accno]=$0
    }
    nseqs++
  } else {
    # the second file

    if (tbd[$1]==1) {
      # This is a sequence which we need.

      # Set a flag that we need to print all lines following this.
      doprint=1
      # We don't need to print this sequence if it comes up again.
      tbd[$1]=0
      # Print the header line
      print header[$1]
      # Done processing this line.
      next
    } else {
      # This is a sequence which we don't need
      doprint=0
    }
  }
}
# If we are currently in the middle of processing a sequence we need, then
# print this line.
doprint{print}

# Summarize results to stderr
END{
  for (i in tbd) ntbd+=tbd[i]
  print ntbd " sequences missing (" nseqs - ntbd -ndups " found) of " nseqs - ndups " total (" ndups " duplicates)" > "/dev/stderr"
}

BEGIN{FS=OFS="\t"} {
  start=$4-1; end=$5; strand=$7;
  size=end-start;
  print $1, start, end, ".", ".", strand, start, end, 0, 1, size, 0
}

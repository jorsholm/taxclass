#!/usr/bin/env awk

BEGIN {
 OFS="\t"
 print "command", "CPU", "time", "mem"
}

/Command being timed:/ {
 sub(/.+Command being timed: /, "")
 command = $0
}

/Percent of CPU this job got:/ {
 sub(/.+Percent of CPU this job got: /, "")
 cpu = $0
}

/Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):/ {
 sub(/.+Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): /, "")
 time = $0
}

/Maximum resident set size \(kbytes\):/ {
 sub(/.+Maximum resident set size \(kbytes\): /, "")
 print command, cpu, time, $0
}

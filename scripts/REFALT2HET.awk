BEGIN {
    FS=" "
    OFS="\t"
}

NR == 1 {
    printf "%s\t%s", $1, $2
    for (i = 3; i <= NF; i += 2) {
        sub(/^REF_/, "", $i)
        printf "\t%s", $i
    }
    printf "\n"
    next
}

{
    printf "%s\t%s", $1, $2
    for (i = 3; i <= NF; i += 2) {
        r = $i + 0  # Force numeric conversion
        a = $(i+1) + 0  # Force numeric conversion
        N = r + a
        if (r == 0 && a == 0) {
            printf "\tNA"
        } else {
            result = ((2 * r * a) / ((r + a) * (r + a))) * (1 + 1/N)
            printf "\t%.6f", result
        }
    }
    printf "\n"
}

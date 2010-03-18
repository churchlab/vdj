{
    if( $0 ~ />.+/ ){
        if( NR > 0 ){ printf "\n" }
        print $1
    } else {
        for( i=1; i<=NF; i++ ){
            printf "%c", $i + 33
        }
    }
}

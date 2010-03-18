{
    if( $0 ~ />.+/ ){
        if( NR > 0 ){ printf "\n" }
        print $1
    } else {
        printf "%s", $1
    }
}

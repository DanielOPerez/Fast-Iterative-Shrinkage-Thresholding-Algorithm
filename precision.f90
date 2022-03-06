  MODULE precision

    ! --------------------------------------------------------------------------
    !  SP: simple precision de la norma IEEE 754
    !  DP: doble precision de la norma IEEE 754
    !
    !  Uso: USE precision, WP => SP o USE precision, WP => DP
    ! --------------------------------------------------------------------------

    INTEGER,PARAMETER::SP=SELECTED_REAL_KIND(6,37)  
    INTEGER,PARAMETER::DP=SELECTED_REAL_KIND(15,307) 

  END MODULE precision

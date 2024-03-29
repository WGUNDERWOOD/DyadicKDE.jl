function kernel(s::Float64, w::Float64, bandwidth::Float64, w_min::Float64,
                w_max::Float64, kernel_name::String)
    if (abs(w - s) <= bandwidth) && (w_min <= w <= w_max) && (w_min <= s <= w_max)
        if kernel_name == "epanechnikov_order_2"
            return kernel_epanechnikov_order_2(s, w, bandwidth, w_min, w_max)

        elseif kernel_name == "epanechnikov_order_4"
            return kernel_epanechnikov_order_4(s, w, bandwidth, w_min, w_max)

        else
            error("Unknown kernel_name")
        end
    end

    return 0.0
end

function kernel_epanechnikov_order_2(s::Float64, w::Float64, bandwidth::Float64,
                                     w_min::Float64, w_max::Float64)
    at = min((w - w_min) / bandwidth, 1)
    bt = min((w_max - w) / bandwidth, 1)
    st = (s - w) / bandwidth

    if (at == 1.0) && (bt == 1.0)
        return 0.75 * (1 / bandwidth) * (1 - st^2)

    else
        c0 = -((48 * (3 * at^4 - 3 * at^3 * bt + at^2 * (-5 + 3 * bt^2) + bt^2 * (-5 + 3 * bt^2) +
                      at * (5 * bt - 3 * bt^3))) /
               ((at + bt)^3 * (60 + 3 * at^4 - 12 * at^3 * bt - 44 * bt^2 + 3 * bt^4 -
                               4 * at * bt * (-8 + 3 * bt^2) + at^2 * (-44 + 30 * bt^2))))
        c1 = -((180 * (at^3 - at^2 * bt + at * (-2 + bt^2) - bt * (-2 + bt^2))) /
               ((at + bt)^3 * (60 + 3 * at^4 - 12 * at^3 * bt - 44 * bt^2 + 3 * bt^4 -
                               4 * at * bt * (-8 + 3 * bt^2) +
                               at^2 * (-44 + 30 * bt^2))))
        c2 = (48 * (3 * at^4 - 3 * at^3 * bt + at^2 * (-5 + 3 * bt^2) + bt^2 * (-5 + 3 * bt^2) +
                    at * (5 * bt - 3 * bt^3))) /
             ((at + bt)^3 * (60 + 3 * at^4 - 12 * at^3 * bt - 44 * bt^2 + 3 * bt^4 -
                             4 * at * bt * (-8 + 3 * bt^2) +
                             at^2 * (-44 + 30 * bt^2)))
        c3 = (180 * (at^3 - at^2 * bt + at * (-2 + bt^2) - bt * (-2 + bt^2))) /
             ((at + bt)^3 * (60 + 3 * at^4 - 12 * at^3 * bt - 44 * bt^2 + 3 * bt^4 -
                             4 * at * bt * (-8 + 3 * bt^2) +
                             at^2 * (-44 + 30 * bt^2)))

        return (c0 + c1 * st + c2 * st^2 + c3 * st^3) / bandwidth
    end
end

function kernel_epanechnikov_order_4(s::Float64, w::Float64, bandwidth::Float64,
                                     w_min::Float64, w_max::Float64)
    at = min((w - w_min) / bandwidth, 1)
    bt = min((w_max - w) / bandwidth, 1)
    st = (s - w) / bandwidth

    if (at == 1.0) && (bt == 1.0)
        return (1 / bandwidth) * ((45 / 32) - (75 / 16) * st^2 + (105 / 32) * st^4)

    else
        c0 = -((240 * (25 * at^12 - 225 * at^11 * bt + 45 * at^10 * (-7 + 25 * bt^2) +
                       at^9 * (2835 * bt - 4125 * bt^3) -
                       30 * at^7 * bt * (252 - 1120 * bt^2 + 705 * bt^4) +
                       15 * at^8 * (56 - 945 * bt^2 + 825 * bt^4) +
                       45 * at^2 * bt^4 * (-588 + 840 * bt^2 - 315 * bt^4 + 25 * bt^6) -
                       9 * at * bt^5 * (-588 + 840 * bt^2 - 315 * bt^4 + 25 * bt^6) +
                       bt^6 * (-588 + 840 * bt^2 - 315 * bt^4 + 25 * bt^6) -
                       15 * at^3 * bt^3 * (-2548 + 4480 * bt^2 - 2240 * bt^4 + 275 * bt^6) +
                       15 * at^4 * bt^2 * (-1764 + 5460 * bt^2 - 4270 * bt^4 + 825 * bt^6) +
                       14 * at^6 * (-42 + 2700 * bt^2 - 4575 * bt^4 + 1775 * bt^6) -
                       6 * at^5 * bt * (-882 + 11200 * bt^2 - 13125 * bt^4 + 3525 * bt^6))) /
               ((at + bt)^7 *
                (8820 + 5 * at^8 - 80 * at^7 * bt - 13160 * bt^2 + 5376 * bt^4 - 540 * bt^6 +
                 5 * bt^8 - 80 * at^5 * bt * (-33 + 26 * bt^2) +
                 20 * at^6 * (-27 + 34 * bt^2) - 16 * at^3 * bt * (651 - 690 * bt^2 + 130 * bt^4) +
                 at^4 * (5376 - 8940 * bt^2 + 3130 * bt^4) -
                 16 * at * bt * (-560 + 651 * bt^2 - 165 * bt^4 + 5 * bt^6) +
                 4 * at^2 * (-3290 + 5334 * bt^2 - 2235 * bt^4 + 170 * bt^6))))
        c1 = -((2100 * (at - bt) *
                (15 * at^10 - 120 * at^9 * bt + at^8 * (-238 + 555 * bt^2) +
                 at^7 * (1904 * bt - 1920 * bt^3) -
                 16 * at^5 * bt * (345 - 679 * bt^2 + 285 * bt^4) +
                 at^6 * (690 - 5656 * bt^2 + 3930 * bt^4) -
                 8 * at * bt^3 * (-504 + 690 * bt^2 - 238 * bt^4 + 15 * bt^6) +
                 bt^4 * (-504 + 690 * bt^2 - 238 * bt^4 + 15 * bt^6) -
                 16 * at^3 * bt * (-252 + 870 * bt^2 - 679 * bt^4 + 120 * bt^6) +
                 at^2 * bt^2 * (-8568 + 13290 * bt^2 - 5656 * bt^4 + 555 * bt^6) +
                 2 * at^4 * (-252 + 6645 * bt^2 - 7798 * bt^4 + 1965 * bt^6))) /
               ((at + bt)^7 *
                (8820 + 5 * at^8 - 80 * at^7 * bt - 13160 * bt^2 + 5376 * bt^4 - 540 * bt^6 +
                 5 * bt^8 - 80 * at^5 * bt * (-33 + 26 * bt^2) +
                 20 * at^6 * (-27 + 34 * bt^2) - 16 * at^3 * bt * (651 - 690 * bt^2 + 130 * bt^4) +
                 at^4 * (5376 - 8940 * bt^2 + 3130 * bt^4) -
                 16 * at * bt * (-560 + 651 * bt^2 - 165 * bt^4 + 5 * bt^6) +
                 4 * at^2 * (-3290 + 5334 * bt^2 - 2235 * bt^4 + 170 * bt^6))))
        c2 = (1200 * (5 * at^12 - 45 * at^11 * bt + 15 * at^10 * (-7 + 15 * bt^2) +
                      at^9 * (945 * bt - 825 * bt^3) -
                      6 * at^7 * bt * (1386 - 2275 * bt^2 + 705 * bt^4) +
                      3 * at^8 * (308 - 1575 * bt^2 + 825 * bt^4) +
                      70 * at^6 * (-35 + 384 * bt^2 - 375 * bt^4 + 71 * bt^6) -
                      30 * at^5 * bt * (-735 + 1778 * bt^2 - 1071 * bt^4 + 141 * bt^6) -
                      9 * at * bt^3 * (1764 - 2450 * bt^2 + 924 * bt^4 - 105 * bt^6 + 5 * bt^8) +
                      bt^4 * (1764 - 2450 * bt^2 + 924 * bt^4 - 105 * bt^6 + 5 * bt^8) +
                      15 * at^2 * bt^2 * (1764 - 3234 * bt^2 + 1792 * bt^4 - 315 * bt^6 + 15 * bt^8) +
                      3 * at^4 * (588 - 16170 * bt^2 + 22680 * bt^4 - 8750 * bt^6 + 825 * bt^8) +
                      at^3 * (-15876 * bt + 59780 * bt^3 - 53340 * bt^5 + 13650 * bt^7 - 825 * bt^9))) /
             ((at + bt)^7 *
              (8820 + 5 * at^8 - 80 * at^7 * bt - 13160 * bt^2 + 5376 * bt^4 - 540 * bt^6 +
               5 * bt^8 -
               80 * at^5 * bt * (-33 + 26 * bt^2) + 20 * at^6 * (-27 + 34 * bt^2) -
               16 * at^3 * bt * (651 - 690 * bt^2 + 130 * bt^4) +
               at^4 * (5376 - 8940 * bt^2 + 3130 * bt^4) -
               16 * at * bt * (-560 + 651 * bt^2 - 165 * bt^4 + 5 * bt^6) +
               4 * at^2 * (-3290 + 5334 * bt^2 - 2235 * bt^4 + 170 * bt^6)))
        c3 = (2100 * (at - bt) *
              (15 * at^10 - 120 * at^9 * bt - 80 * at^7 * bt * (-25 + 24 * bt^2) +
               5 * at^8 * (-50 + 111 * bt^2) - 16 * at^5 * bt * (462 - 775 * bt^2 + 285 * bt^4) +
               at^6 * (924 - 6100 * bt^2 + 3930 * bt^4) -
               16 * at^3 * bt * (-630 + 1302 * bt^2 - 775 * bt^4 + 120 * bt^6) +
               2 * at^4 * (-630 + 8274 * bt^2 - 8650 * bt^4 + 1965 * bt^6) -
               8 * at * bt * (588 - 1260 * bt^2 + 924 * bt^4 -
                              250 * bt^6 + 15 * bt^8) +
               bt^2 * (588 - 1260 * bt^2 + 924 * bt^4 - 250 * bt^6 + 15 * bt^8) +
               at^2 * (588 - 12600 * bt^2 + 16548 * bt^4 - 6100 * bt^6 + 555 * bt^8))) /
             ((at + bt)^7 *
              (8820 + 5 * at^8 - 80 * at^7 * bt - 13160 * bt^2 + 5376 * bt^4 - 540 * bt^6 +
               5 * bt^8 - 80 * at^5 * bt * (-33 + 26 * bt^2) +
               20 * at^6 * (-27 + 34 * bt^2) - 16 * at^3 * bt * (651 - 690 * bt^2 + 130 * bt^4) +
               at^4 * (5376 - 8940 * bt^2 + 3130 * bt^4) -
               16 * at * bt * (-560 + 651 * bt^2 - 165 * bt^4 + 5 * bt^6) +
               4 * at^2 * (-3290 + 5334 * bt^2 - 2235 * bt^4 + 170 * bt^6)))
        c4 = (3360 * (15 * at^10 - 135 * at^9 * bt + 135 * at^8 * (-2 + 5 * bt^2) -
                      45 * at^7 * bt * (-54 + 55 * bt^2) -
                      3 * at^5 * bt * (2499 - 4750 * bt^2 +
                                       1950 * bt^4) + at^6 * (833 - 6900 * bt^2 + 4800 * bt^4) -
                      9 * at * bt^3 * (-630 + 833 * bt^2 - 270 * bt^4 + 15 * bt^6) +
                      bt^4 * (-630 + 833 * bt^2 - 270 * bt^4 + 15 * bt^6) +
                      15 * at^2 * bt^2 * (-630 + 1029 * bt^2 - 460 * bt^4 + 45 * bt^6) +
                      15 * at^4 * (-42 + 1029 * bt^2 - 1230 * bt^4 + 320 * bt^6) -
                      5 * at^3 * bt * (-1134 + 3724 * bt^2 - 2850 * bt^4 + 495 * bt^6))) /
             ((at + bt)^7 *
              (8820 + 5 * at^8 - 80 * at^7 * bt - 13160 * bt^2 + 5376 * bt^4 - 540 * bt^6 +
               5 * bt^8 - 80 * at^5 * bt * (-33 + 26 * bt^2) +
               20 * at^6 * (-27 + 34 * bt^2) - 16 * at^3 * bt * (651 - 690 * bt^2 + 130 * bt^4) +
               at^4 * (5376 - 8940 * bt^2 + 3130 * bt^4) -
               16 * at * bt * (-560 + 651 * bt^2 - 165 * bt^4 + 5 * bt^6) +
               4 * at^2 * (-3290 + 5334 * bt^2 - 2235 * bt^4 + 170 * bt^6)))
        c5 = (12600 * (at - bt) *
              (2 * at^8 - 16 * at^7 * bt - 8 * at^5 * bt * (-39 + 32 * bt^2) +
               at^6 * (-39 + 74 * bt^2) -
               16 * at^3 * bt * (63 - 72 * bt^2 + 16 * bt^4) +
               at^4 * (126 - 543 * bt^2 + 284 * bt^4) -
               8 * at * bt * (-98 + 126 * bt^2 - 39 * bt^4 + 2 * bt^6) +
               bt^2 * (-98 + 126 * bt^2 - 39 * bt^4 + 2 * bt^6) +
               at^2 * (-98 + 672 * bt^2 - 543 * bt^4 + 74 * bt^6))) /
             ((at + bt)^7 *
              (8820 + 5 * at^8 - 80 * at^7 * bt - 13160 * bt^2 + 5376 * bt^4 - 540 * bt^6 +
               5 * bt^8 - 80 * at^5 * bt * (-33 + 26 * bt^2) +
               20 * at^6 * (-27 + 34 * bt^2) - 16 * at^3 * bt * (651 - 690 * bt^2 + 130 * bt^4) +
               at^4 * (5376 - 8940 * bt^2 + 3130 * bt^4) -
               16 * at * bt * (-560 + 651 * bt^2 - 165 * bt^4 + 5 * bt^6) +
               4 * at^2 * (-3290 + 5334 * bt^2 - 2235 * bt^4 + 170 * bt^6)))

        return (c0 + c1 * st + c2 * st^2 + c3 * st^3 + c4 * st^4 + c5 * st^5) / bandwidth
    end
end

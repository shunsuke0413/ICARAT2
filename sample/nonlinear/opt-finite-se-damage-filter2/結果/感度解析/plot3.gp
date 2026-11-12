# Sensitivity comparison (gnuplot)
# Run: gnuplot plot3.gp (or load "plot3.gp" in gnuplot).
# Batch/PDF only: gnuplot -e "show_plots=0" plot3.gp
if (!exists("show_plots")) show_plots = 1

# Legend placement. All panels default to bottom right on each load.
# Edit the settings below to change positions.
# Examples: "top left", "bottom right", "outside top center horizontal"
# Fine placement: "at graph 0.98,0.98 right top" (graph coordinates: 0 to 1).
# Example setting below: legend_position = "bottom left"
# Each panel can also be positioned separately via the four variables below.
legend_position = "bottom right"
legend_sensitivity = ""
legend_error = ""
legend_agreement = ""
legend_volume = ""
# Empty panel settings follow legend_position, including subsequent reloads.
legend(panel) = (strlen(panel) > 0 ? panel : legend_position)
unset logscale
set autoscale
unset multiplot
unset output
file = "/Users/tsurugashunsuke/lab/ICARAT2/sample/nonlinear/opt-finite-se-damage-filter2/sensitivity_comparison.csv"  # CSV path
set datafile separator comma
set terminal pdfcairo enhanced font "Helvetica,14" size 16cm,10cm
set grid
eval "set key " . legend_position
blue = "#0072BD"
orange = "#D95319"
red = "#C00000"

max(a,b) = (a > b ? a : b)

# Select rows by perturbation size; NaN means 'do not plot'.
match(h) = abs(column("h")-h) < h*0.001
ok(h) = match(h) && strcol("branch_status") ne "BRANCH-CHANGE"
bad(h) = match(h) && strcol("branch_status") eq "BRANCH-CHANGE"
# 1, 3, 4: one figure for each perturbation size
 do for [i=1:2] {
    h = (i==1 ? 1e-5 : 1e-6)
    tag = (i==1 ? "1e-5" : "1e-6")
    # Filter rows before drawing lines: NaN on every other input row
    # otherwise breaks the curve into isolated, invisible line vertices.
    set table $curve separator comma
    plot file using (column("element")):(sprintf("%.17g",column("adjoint"))):(sprintf("%.17g",column("volume_adjoint"))):(sprintf("%.17g",column("finite_difference"))):(sprintf("%.17g",column("volume_finite_difference"))) with table if (match(h))
    unset table
    set xlabel "Element"
    set ylabel "Objective sensitivity"
    eval "set key " . legend(legend_sensitivity)
    set output sprintf("sensitivity_%s.pdf",tag)
    plot $curve using 1:2 with lines dt 1 lw 2 lc rgb blue title "Adjoint", \
         $curve using 1:4 with lines dt 2 lw 2 lc rgb orange title "Finite difference", \
         file using (bad(h)?column("element"):1/0):(column("finite_difference")) with points pt 2 ps 1.5 lw 2 lc rgb red notitle

    set xlabel "Adjoint sensitivity"
    set ylabel "Finite difference sensitivity"
    eval "set key " . legend(legend_agreement)
    set output sprintf("agreement_%s.pdf",tag)
    plot x with lines dt 2 lc rgb "#777777" title "y=x", \
         file using (ok(h)?column("adjoint"):1/0):(column("finite_difference")) with points pt 7 ps 0.8 lc rgb blue title "Smooth branches", \
         file using (bad(h)?column("adjoint"):1/0):(column("finite_difference")) with points pt 2 ps 1.5 lw 2 lc rgb red notitle

    set xlabel "Element"
    set ylabel "Volume constraint sensitivity"
    eval "set key " . legend(legend_volume)
    set output sprintf("volume_%s.pdf",tag)
    plot $curve using 1:3 with lines dt 1 lw 2 lc rgb blue title "Analytical", \
         $curve using 1:5 with lines dt 2 lw 2 lc rgb orange title "Finite difference", \
         file using (bad(h)?column("element"):1/0):(column("volume_finite_difference")) with points pt 2 ps 1.5 lw 2 lc rgb red notitle
}

# 2: relative error, logarithmic y-axis (zero shown at 1e-12)
set xlabel "Element"
set ylabel "Relative error"
eval "set key " . legend(legend_error)
set logscale y
set yrange [1e-12:*]
set output "relative_error.pdf"
plot file using (ok(1e-5)?column("element"):1/0):(max(1e-12,column("relative_error"))) with points pt 7 lc rgb blue title "h=1e-5", \
     file using (ok(1e-6)?column("element"):1/0):(max(1e-12,column("relative_error"))) with points pt 5 lc rgb orange title "h=1e-6", \
     file using (bad(1e-5)?column("element"):1/0):(max(1e-12,column("relative_error"))) with points pt 2 ps 1.5 lc rgb red notitle, \
     file using (bad(1e-6)?column("element"):1/0):(max(1e-12,column("relative_error"))) with points pt 4 ps 1.5 lc rgb red notitle
unset output

unset logscale y
set autoscale

# Display the CSV comparison as well as writing the PDFs above.
# qt persist keeps each window visible after this script has finished.
if (show_plots) {
    do for [i=1:2] {
        h = (i==1 ? 1e-5 : 1e-6)
        tag = (i==1 ? "1e-5" : "1e-6")
    # Filter rows before drawing lines: NaN on every other input row
    # otherwise breaks the curve into isolated, invisible line vertices.
    set table $curve separator comma
    plot file using (column("element")):(sprintf("%.17g",column("adjoint"))):(sprintf("%.17g",column("volume_adjoint"))):(sprintf("%.17g",column("finite_difference"))):(sprintf("%.17g",column("volume_finite_difference"))) with table if (match(h))
    unset table
        set terminal qt (i-1) persist enhanced font "Helvetica,11" size 1200,800 title sprintf("Sensitivity comparison: h=%s",tag)
        unset output
        set multiplot layout 2,2 title sprintf("Adjoint versus finite differences (h=%s)",tag)
        set autoscale
        unset logscale
        eval "set key " . legend_position
        set xlabel "Element"
        set ylabel "Objective sensitivity"
        eval "set key " . legend(legend_sensitivity)
        set title "Sensitivity by element"
        plot $curve using 1:2 with lines dt 1 lw 2 lc rgb blue title "Adjoint", \
             $curve using 1:4 with lines dt 2 lw 2 lc rgb orange title "Finite difference", \
             file using (bad(h)?column("element"):1/0):(column("finite_difference")) with points pt 2 ps 1.2 lc rgb red notitle

        set title "Relative error (fraction)"
        set ylabel "Relative error"
        eval "set key " . legend(legend_error)
        set logscale y
        set yrange [1e-12:*]
        plot file using (ok(h)?column("element"):1/0):(max(1e-12,column("relative_error"))) with points pt 7 ps 0.6 lc rgb blue title "Smooth branches", \
             file using (bad(h)?column("element"):1/0):(max(1e-12,column("relative_error"))) with points pt 2 ps 1.2 lc rgb red notitle

        unset logscale y
        set autoscale
        # Explicit shared limits prevent the y=x function inheriting element indices.
        stats file using (match(h)?column("adjoint"):1/0):(column("finite_difference")) nooutput
        lower = (STATS_min_x < STATS_min_y ? STATS_min_x : STATS_min_y)
        upper = (STATS_max_x > STATS_max_y ? STATS_max_x : STATS_max_y)
        margin = max(1e-12,0.05*(upper-lower))
        set xrange [lower-margin:upper+margin]
        set yrange [lower-margin:upper+margin]
        set title "Agreement (ideal: y=x)"
        set xlabel "Adjoint sensitivity"
        set ylabel "Finite difference sensitivity"
        eval "set key " . legend(legend_agreement)
        plot x with lines dt 2 lc rgb "#777777" title "y=x", \
             file using (ok(h)?column("adjoint"):1/0):(column("finite_difference")) with points pt 7 ps 0.6 lc rgb blue title "Smooth branches", \
             file using (bad(h)?column("adjoint"):1/0):(column("finite_difference")) with points pt 2 ps 1.2 lc rgb red notitle

        set autoscale
        set title "Volume constraint sensitivity"
        set xlabel "Element"
        set ylabel "Volume constraint sensitivity"
        eval "set key " . legend(legend_volume)
        plot $curve using 1:3 with lines dt 1 lw 2 lc rgb blue title "Analytical", \
             $curve using 1:5 with lines dt 2 lw 2 lc rgb orange title "Finite difference", \
             file using (bad(h)?column("element"):1/0):(column("volume_finite_difference")) with points pt 2 ps 1.2 lc rgb red notitle
        unset multiplot
    }
}

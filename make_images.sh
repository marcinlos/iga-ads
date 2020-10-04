#!/bin/bash


for i in `seq 0 100 3000`; do
    # tumor
    input_file="data/tumor_${i}.data"
    output_file="img/tumor_${i}.png"
    value_range="1.5"
    title="Cancer_cell_density"

    gnuplot -e "input_file='"${input_file}"';\
                output_file='"${output_file}"';\
                value_range="${value_range}";\
                plot_title='"${title}"'" plot

    input_file="data/taf_${i}.data"
    output_file="img/taf_${i}.png"
    value_range="1"
    title="TAF"

    gnuplot -e "input_file='"${input_file}"';\
                output_file='"${output_file}"';\
                value_range="${value_range}";\
                plot_title='"${title}"'" plot

    #input_file="data/ecm_${i}.data"
    #output_file="img/ecm_${i}.png"
    #value_range="1"
    #title="ECM"

    #gnuplot -e "input_file='"${input_file}"';\
    #            output_file='"${output_file}"';\
    #            value_range="${value_range}";\
    #            plot_title='"${title}"'" plot

    #input_file="data/degraded_ecm_${i}.data"
    #output_file="img/degraded_ecm_${i}.png"
    #value_range="0.04"
    #title="Degraded_ECM"

    #gnuplot -e "input_file='"${input_file}"';\
    #            output_file='"${output_file}"';\
    #            value_range="${value_range}";\
    #            plot_title='"${title}"'" plot

    #input_file="data/fibronectin_${i}.data"
    #output_file="img/fibronectin_${i}.png"
    #value_range="1"
    #title="Fibronectin"

    #gnuplot -e "input_file='"${input_file}"';\
    #            output_file='"${output_file}"';\
    #            value_range="${value_range}";\
    #            plot_title='"${title}"'" plot

    input_file="data/vasculature_${i}.data"
    output_file="img/vasculature_${i}.png"
    value_range="1"
    title="Vasculature"

    gnuplot -e "input_file='"${input_file}"';\
                output_file='"${output_file}"';\
                value_range="${value_range}";\
                plot_title='"${title}"'" plot

    input_file="data/oxygen_${i}.data"
    output_file="img/oxygen_${i}.png"
    value_range="1"
    title="Oxygen"

    gnuplot -e "input_file='"${input_file}"';\
                output_file='"${output_file}"';\
                value_range="${value_range}";\
                plot_title='"${title}"'" plot
done

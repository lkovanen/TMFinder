## Plot the temporal motifs in the small test data.

## Description of node types, including node colors.
type_desc="type_desc_small.txt"

plot_code="../python/plot_motifs.py"

nodes="42 66"
events="1"

N_v=4
N_h=8
order_by='count'
group_by="untyped motif"
label_def='{count:d}'
min_count=10 # Motif must have a count of at least min_count to be included.

test_output="test_small_output.dat"

# Plot all motifs ordered by motif count.
python ${plot_code} --most_common -df ${test_output} -tf ${type_desc} --prefix ${test_output/.dat} -nt ${nodes} -et ${events} --verbose --label ${label_def} --order_by ${order_by} -or

# Plot all 2-event motifs, each topology in a separate figure
python ${plot_code} --most_common_grouped_single -e 2 -df ${test_output} -tf ${type_desc} --prefix ${test_output/.dat} -nt ${nodes} -et ${events} --verbose --label ${label_def} --group_by ${group_by} --order_by ${order_by} -or

# Plot all motifs in a single figure, grouped by topology and sorted
# by count in each group.
python ${plot_code} --most_common_grouped -df ${test_output} -tf ${type_desc} --prefix ${test_output/.dat} -nt ${nodes} -et ${events} --verbose --label ${label_def} --group_by ${group_by} --order_by ${order_by} -or

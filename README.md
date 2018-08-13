# LDzoom
An R script to reduce a hapmap based on LD decay around a specified marker

Usage: LDzoom.R [options] > marker.positions

Options:
	--marker=CHARACTER
		Marker to be used a center of the interval and to compute LD around it.

	--hmp=CHARACTER
		Hapmap to compute the LD and zoom into it.

	--chrom=CHARACTER
		Chromosome where --marker is located (both --chrom, --start and --end are used to specify which region to zoom in for those markers with multiple mappings to different loci.

	--start=CHARACTER
		Start position of --marker alignment to the reference (both --chrom, --start and --end are used to specify which region to zoom in for those markers with multiple mappings to different loci.

	--end=CHARACTER
		End position of --marker alignment to the reference (both --chrom, --start and --end are used to specify which region to zoom in for those markers with multiple mappings to different loci.

	--ldthres=CHARACTER
		When the LD from a given position to --marker goes below this --ldthres, the previous position is a end limit of the final interval.

	--interval_max=CHARACTER
		The final interval reported will never exceed --interval_max, even if the LD does not ever go below --ldthres

	--interval_min=CHARACTER
		The final interval span will never be smaller than --interval_min, even if the LD goes below --ldthres closer to the --marker

	--window_size=CHARACTER
		How many adjacent markers will be used to compute each LD value.

	--window_step=CHARACTER
		How many markers will skip one LD window from the previous LD window.

	-h, --help
		Show this help message and exit

Note: marker, hmp, chrom, start and end are mandatory parameters.
marker, chrom, start and end define a marker locus, which must be included in the region defined by markers in hmp.
ldthres value will help define the new region based on LD decay.

Output: the script generates several files, within a LD/ folder which must exist.
Also outputs to the stdout the a row with the following format:

#marker_id:marker_chrom-marker_start:marker_LD_start-marker_LD_end

marker_id marker_chrom and marker_start are the same used as parameters
marker_LD_start and marker_LD_end are the positions defined by LD decay and the specified --ldthress

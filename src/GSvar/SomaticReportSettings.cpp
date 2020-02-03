#include "SomaticReportSettings.h"
#include "SomaticReportHelper.h"


SomaticReportSettings::SomaticReportSettings()
	: report_config()
	, tumor_ps()
	, normal_ps()
	, filters()
	, report_type()
	, include_gap_statistics(false)
{
}

VariantList SomaticReportSettings::filterVariants(const VariantList &snvs, const SomaticReportSettings& sett)
{
	QSet<int> variant_indices = sett.report_config.variantIndices(VariantType::SNVS_INDELS, false).toSet();

	VariantList result;

	result.copyMetaData(snvs);

	FilterResult filter_res =sett.filters.apply(snvs); //does not regard "include" result of report_config

	//Adapt filter results to results from report settings
	for(int index : variant_indices)
	{
		filter_res.flags()[index] = sett.report_config.variantConfig(index).showInReport();
	}


	result.addAnnotation("alt_var_alteration","If an alternative text for protein change is specified in report config, this is stored here.", "");
	result.addAnnotation("alt_var_description", "Alternate description text for variant alteration", "");
	for(int i=0; i<snvs.count(); ++i)
	{
		if(!filter_res.flags()[i]) continue;

		result.append(snvs[i]);

		//add additional report config info into new empty annotation columns
		if(variant_indices.contains(i) && sett.report_config.variantConfig(i).showInReport())
		{
			result[result.count()-1].annotations().append(sett.report_config.variantConfig(i).include_variant_alteration.toUtf8());
			result[result.count()-1].annotations().append(sett.report_config.variantConfig(i).include_variant_description.toUtf8());
		}
		else
		{
			result[result.count()-1].annotations().append({"",""}); //empty annotation columns
		}
	}

	return result;
}

CnvList SomaticReportSettings::filterCnvs(const CnvList &cnvs, const SomaticReportSettings &sett)
{
	QSet<int> cnv_indices = sett.report_config.variantIndices(VariantType::CNVS, false).toSet();

	CnvList result;
	result.copyMetaData(cnvs);

	QBitArray cnv_flags(cnvs.count(), true);

	for(int index : cnv_indices)
	{
		cnv_flags[index] = sett.report_config.variantConfig(index).showInReport();
	}

	for(int i=0; i<cnvs.count(); ++i)
	{
		if(!cnv_flags[i]) continue;

		result.append(cnvs[i]);
	}
	return result;
}

#include "SomaticReportConfiguration.h"

SomaticReportVariantConfiguration::SomaticReportVariantConfiguration()
	: variant_type(VariantType::SNVS_INDELS)
	, variant_index(-1)
	, report_type()
	, exclude_artefact(false)
	, exclude_low_tumor_content(false)
	, exclude_low_copy_number(false)
	, exclude_high_baf_deviation(false)
	, exclude_other_reason(false)
	, include_variant_alteration()
	, include_variant_description()
	, comment()
{
}

bool SomaticReportVariantConfiguration::showInReport() const
{
	return !(exclude_artefact || exclude_low_tumor_content || exclude_low_copy_number || exclude_high_baf_deviation || exclude_other_reason);
}

SomaticReportConfiguration::SomaticReportConfiguration()
	: variant_config_()
	, created_by_(Helper::userName())
	, created_at_(QDateTime::currentDateTime())
{
}

const QList<SomaticReportVariantConfiguration>& SomaticReportConfiguration::variantConfig() const
{
	return variant_config_;
}

QList<int> SomaticReportConfiguration::variantIndices(VariantType type, bool only_selected, QString report_type) const
{
	QList<int> output;

	for(const auto& var_conf : variant_config_)
	{
		if (var_conf.variant_type!=type) continue;
		if (only_selected && !var_conf.showInReport()) continue;
		if (!report_type.isNull() && var_conf.report_type!=report_type) continue;
		output << var_conf.variant_index;
	}

	std::sort(output.begin(), output.end());

	return output;
}


bool SomaticReportConfiguration::exists(VariantType type, int index) const
{
	for(const auto& var_conf : variant_config_)
	{
		if (var_conf.variant_index==index && var_conf.variant_type==type) return true;
	}

	return false;
}

bool SomaticReportConfiguration::set(const SomaticReportVariantConfiguration &config)
{
	//set variant config (if already contained in list)
	for (int i=0; i<variant_config_.count(); ++i)
	{
		const SomaticReportVariantConfiguration& var_conf = variant_config_[i];
		if (var_conf.variant_index==config.variant_index && var_conf.variant_type==config.variant_type)
		{
			variant_config_[i] = config;
			return true;
		}
	}
	//set variant config (if not yet contained)
	variant_config_ << config;

	sortByPosition();

	return false;
}

bool SomaticReportConfiguration::remove(VariantType type, int index)
{
	for(int i=0; i<variant_config_.count(); ++i)
	{
		const auto& var_conf = variant_config_[i];
		if (var_conf.variant_index==index && var_conf.variant_type==type)
		{
			variant_config_.removeAt(i);
			return true;
		}
	}
	return false;
}

int SomaticReportConfiguration::count()
{
	return variant_config_.count();
}

QDateTime SomaticReportConfiguration::createdAt() const
{
	return created_at_;
}
QString SomaticReportConfiguration::createdBy() const
{
	return created_by_;
}
void SomaticReportConfiguration::setCreatedAt(QDateTime time)
{
	created_at_ = time;
}
void SomaticReportConfiguration::setCreatedBy(QString user)
{
	created_by_ = user;
}
void SomaticReportConfiguration::sortByPosition()
{
	std::sort(variant_config_.begin(), variant_config_.end(), [](const SomaticReportVariantConfiguration& a, const SomaticReportVariantConfiguration& b){return a.variant_index < b.variant_index;});
}

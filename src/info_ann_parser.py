def parse_info_ann(info_ann_str: str, ann_field_info: dict) -> list:
    annotations = []
    ann_description: str = ann_field_info.get('Description', None)
    if ann_description and info_ann_str:
        ann_keys = ann_description.split(": ")[1].replace("'","").split("|")
        for ann_str in info_ann_str.split(','):
            ann_values = ann_str.split('|')
            ann_dict = dict()
            for i in range(len(ann_keys)):
                ann_key = ann_keys[i].strip()
                ann_value = ann_values[i].strip()
                ann_dict[ann_key] = ann_value
            annotations.append(ann_dict)

    return annotations
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

def remove_fields_from_bib(bib_file_path, output_path, fields_to_remove):
    with open(bib_file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    for field in fields_to_remove:
        pattern = re.compile(r'\s*'+ re.escape(field) + r'\s*=\s*\{.*?\}\s*,?', re.DOTALL)
        content = re.sub(pattern, '', content)

    with open(output_path, 'w', encoding='utf-8') as file:
        file.write(content)

bib_file_path = 'references.bib'
output_path = 'ref.bib'
fields_to_remove = ['note', 'url', 'urldate']
remove_fields_from_bib(bib_file_path, output_path, fields_to_remove)

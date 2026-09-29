#!/usr/bin/env python
# -*- coding: utf-8 -*-
from docx import Document

doc = Document(r'D:\IA_AAA_analysis\Manuscript_v4.docx')

print('=== FULL DOCUMENT DUMP ===')
for i, para in enumerate(doc.paragraphs):
    text = para.text.strip()
    if text:
        style = para.style.name
        print(f'[{i}] ({style}) {text[:200]}')
        if len(text) > 200:
            print(f'    ... {text[200:400]}')

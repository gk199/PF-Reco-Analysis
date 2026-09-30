"""Select inputs and write a deterministic per-file job table."""
import hashlib
import os
from pathlib import Path, PurePosixPath
import shlex
import sys


def select(lines, host, offset, count):
    urls, seen = [], set()
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if any(c.isspace() for c in line):
            raise ValueError('Unexpected whitespace in input URL: ' + line)
        if line.startswith('/store/'):
            line = host.rstrip('/') + '/' + line
        if not line.startswith('root://'):
            raise ValueError('Expected /store/... or root://...: ' + line)
        if line not in seen:
            seen.add(line)
            urls.append(line)
    urls = urls[offset:]
    if count:
        urls = urls[:count]
    if not urls:
        raise ValueError('No inputs selected; check INPUT_LIST / INPUT_FILES and FILE_OFFSET.')
    return urls


if __name__ == '__main__':
    source, count, offset, host, output = sys.argv[1:]
    explicit = os.environ.get('INPUT_FILES', '').strip()
    lines = shlex.split(explicit) if explicit else Path(source).read_text().splitlines()
    urls = select(lines, host, int(offset), int(count))
    out = Path(output)
    (out / 'selected_inputs.txt').write_text('\n'.join(urls) + '\n')
    rows = []
    for i, url in enumerate(urls):
        tag = f'{i:05d}_{PurePosixPath(url).stem}_{hashlib.sha256(url.encode()).hexdigest()[:10]}'
        if not all(c.isalnum() or c in '_-.' for c in tag):
            raise ValueError('Unsupported characters in input basename: ' + tag)
        rows.append(f'{i} {tag}')
    (out / 'jobs.txt').write_text('\n'.join(rows) + '\n')
    print(f'Selected {len(urls)} file(s), starting at offset {offset}.')
    for url in urls[:5]:
        print(url)

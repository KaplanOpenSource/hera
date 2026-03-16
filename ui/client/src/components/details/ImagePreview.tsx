import { Box } from '@mui/material';
import { useEffect, useState } from 'react';
import { fetchPython } from '../../io/fetchPython';

const IMAGE_EXTENSIONS = /\.(png|jpg|jpeg|bmp|gif|tiff|tif|webp)$/i;

export const isImageResource = (resource: unknown): resource is string => {
  return typeof resource === 'string' && IMAGE_EXTENSIONS.test(resource);
};

export const ImagePreview = ({
  path,
}: {
  path: string;
}) => {
  const [dataUrl, setDataUrl] = useState<string | null>(null);

  useEffect(() => {
    setDataUrl(null);
    const ext = path.match(IMAGE_EXTENSIONS)?.[1]?.toLowerCase() ?? 'png';
    const mime = ext === 'jpg' ? 'image/jpeg'
      : ext === 'tif' || ext === 'tiff' ? 'image/tiff'
      : `image/${ext}`;

    fetchPython({
      results: ['b64'],
      code: `
import base64
with open('${path}', 'rb') as f:
    b64 = base64.b64encode(f.read()).decode('ascii')
`,
    }).then(({ data, problem }) => {
      if (!problem && data?.b64) {
        setDataUrl(`data:${mime};base64,${data.b64}`);
      }
    });
  }, [path]);

  return dataUrl
    ? (
      <Box sx={{ height: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', overflow: 'auto' }}>
        <img
          src={dataUrl}
          style={{ maxWidth: '100%', maxHeight: '100%', objectFit: 'contain' }}
        />
      </Box>
    )
    : null;
};

import { Box } from '@mui/material';
import { BASEURL } from '../../shared/baseurl';

const IMAGE_EXTENSIONS = /\.(png|jpg|jpeg|bmp|gif|tiff|tif|webp)$/i;

export const isImageResource = (resource: unknown): resource is string => {
  return typeof resource === 'string' && IMAGE_EXTENSIONS.test(resource);
};

export const ImagePreview = ({
  path,
}: {
  path: string;
}) => {
  const url = `${BASEURL}/file${path.startsWith('/') ? '' : '/'}${path}`;

  return (
    <Box sx={{ height: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', overflow: 'auto' }}>
      <img
        src={url}
        style={{ maxWidth: '100%', maxHeight: '100%', objectFit: 'contain' }}
      />
    </Box>
  );
};

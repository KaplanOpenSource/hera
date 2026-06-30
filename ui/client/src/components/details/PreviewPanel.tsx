import { useEffect, useState } from 'react';
import { fetchDocument } from '../../io/FetchDocument';
import { ProjectDocument } from '../../shared/types';
import { ImagePreview, isImageResource } from './ImagePreview';
import { GeoJsonPreview, isGeoJsonResource } from './maps/GeoJsonPreview';
import { isTileUrl, TileMapView } from './maps/TileMapView';

export const hasPreview = (doc: ProjectDocument): boolean => {
  return isTileUrl(doc.resource) || isImageResource(doc.resource) || isGeoJsonResource(doc.resource);
};

export const PreviewPanel = ({
  docid,
}: {
  docid: string;
}) => {
  const [resource, setResource] = useState<string | undefined>(undefined);

  useEffect(() => {
    (async () => {
      if (docid) {
        const data = await fetchDocument(docid);
        setResource(data?.resource);
      } else {
        setResource(undefined);
      }
    })();
  }, [docid]);

  if (!resource) {
    return null;
  }
  if (isTileUrl(resource)) {
    return <TileMapView url={resource} />;
  }
  if (isGeoJsonResource(resource)) {
    return <GeoJsonPreview path={resource} />;
  }
  if (isImageResource(resource)) {
    return <ImagePreview path={resource} />;
  }
  return null;
};

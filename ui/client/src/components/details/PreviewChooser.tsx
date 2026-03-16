import { ProjectDocument } from '../../shared/types';
import { ImagePreview, isImageResource } from './ImagePreview';
import { GeoJsonPreview, isGeoJsonResource } from './maps/GeoJsonPreview';
import { isTileUrl, TileMapView } from './maps/TileMapView';

export const hasPreview = (doc: ProjectDocument): boolean => {
  return isTileUrl(doc.resource) || isImageResource(doc.resource) || isGeoJsonResource(doc.resource);
};

export const PreviewChooser = ({
  doc,
}: {
  doc: ProjectDocument;
}) => {
  if (isTileUrl(doc.resource)) {
    return <TileMapView url={doc.resource} />;
  }
  if (isGeoJsonResource(doc.resource)) {
    return <GeoJsonPreview path={doc.resource} />;
  }
  if (isImageResource(doc.resource)) {
    return <ImagePreview path={doc.resource} />;
  }
  return null;
};

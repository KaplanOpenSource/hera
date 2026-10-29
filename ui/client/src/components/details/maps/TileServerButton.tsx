import { useState } from 'react';
import { Autocomplete, Box, Button, DialogActions, DialogContent, DialogTitle, TextField, Typography } from '@mui/material';
import { Layers } from '@mui/icons-material';
import { documentApi } from '../../../api/api';
import { ButtonDialog } from '../../../elements/ButtonDialog';
import { Collection } from '../../../shared/types';
import { useProjectStore } from '../../../stores/useProjectStore';
import { findConfigDoc } from './projectConfigDoc';
import { TILE_SERVERS } from './tileServers';

export const TileServerButton = ({
  currentUrl,
}: {
  currentUrl: string,
}) => {
  const [url, setUrl] = useState(currentUrl);

  const save = async () => {
    const project = useProjectStore.getState().currProject;
    if (!project) {
      return;
    }
    const config = findConfigDoc(project);
    if (config) {
      const saved = await documentApi.update({ ...config, desc: { ...config.desc, defaultTileServer: url } });
      if (saved) {
        useProjectStore.getState().setCurrentProject(saved.project);
      }
      return;
    }
    const created = await documentApi.create({
      projectName: project.name,
      collection: Collection.Cache,
      type: `${project.name}__config__`,
      resource: '',
      dataFormat: 'string',
      desc: { defaultTileServer: url },
    }, 'set default tile server');
    if (created) {
      useProjectStore.getState().setCurrentProject(created.project);
    }
  };

  return (
    <ButtonDialog
      icon={<Layers fontSize="small" />}
      title="Default tile server"
      onOpen={() => setUrl(currentUrl)}
      sx={{ bgcolor: 'primary.main', color: 'primary.contrastText', boxShadow: 2, '&:hover': { bgcolor: 'primary.dark' } }}
      dialogProps={{ fullWidth: true, maxWidth: 'sm' }}
    >
      {(close) => (
        <>
          <DialogTitle>Default tile server</DialogTitle>
          <DialogContent>
            <Typography sx={{ mb: 1 }}>
              Pick a tile server from the list, or type the URL of your own tile server.
              It is saved as defaultTileServer in the project config.
            </Typography>
            <Autocomplete
              freeSolo
              disableClearable
              forcePopupIcon
              openOnFocus
              options={TILE_SERVERS.map(server => server.url)}
              filterOptions={options => options}
              inputValue={url}
              onInputChange={(_e, value) => setUrl(value)}
              renderOption={(props, option) => (
                <Box component="li" {...props} key={option} sx={{ flexDirection: 'column', alignItems: 'flex-start !important' }}>
                  <Typography variant="body2">{TILE_SERVERS.find(server => server.url === option)?.name}</Typography>
                  <Typography variant="caption" color="text.secondary">{option}</Typography>
                </Box>
              )}
              renderInput={params => (
                <TextField
                  {...params}
                  autoFocus
                  margin="dense"
                  size="small"
                  label="Tile URL"
                />
              )}
            />
          </DialogContent>
          <DialogActions>
            <Button onClick={close}>Cancel</Button>
            <Button
              variant="contained"
              disabled={!url.trim()}
              onClick={async () => {
                await save();
                close();
              }}
            >
              Save
            </Button>
          </DialogActions>
        </>
      )}
    </ButtonDialog>
  );
};

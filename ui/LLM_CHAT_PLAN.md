# LLM chat panel - plan

Issue: #543 (UI - space for internal LLM)

## Goal

A chat panel in the UI. The user asks about the open project.
Example: "This is a risk project, what do you suggest to do now?"
The LLM answers, and can suggest Hera documents to create.
The user approves each document before it is created.

## Main idea

- The LLM only writes text. It never touches the database.
- The Hera server adds the project details to each question.
- The Hera server turns suggestions into real documents, only after the user clicks.

Flow:

```
UI chat -> Hera server (adds project details) -> RAG (docs + code) -> Ollama (LLM)
                                                                          |
UI shows answer + suggestion cards <- Hera server (checks the JSON) <-----+
UI "Create" click -> Hera server -> normal Hera code (addCacheDocument, ...)
```

## Step 0 - Prep: run the LLM + RAG with our stack

The RAG code is in `hera/utils/rag/`. Ollama is already in `docker-compose.yml`.

1. Add Qdrant, Cassandra and the RAG server to `docker-compose.yml`, under the `ai` profile.
2. Build the RAG index (docs + code).
3. Ask it real project questions. Check the answers.
4. If llama3 is weak, try a better model before going on.

Done when: `docker compose --profile ai up` starts everything, and answers are useful.

## Step 1 - Chat panel, answers only

1. Server: a chat endpoint. It collects the project details (toolkits, documents, settings).
   It sends them with the question to the RAG, and streams the answer back.
2. RAG: accept extra context (the project details) next to the question.
3. Client: a chat panel on the right side. It shows the answer as it streams.

Done when: the user can ask about the project and get a relevant answer.

## Step 2 - Chat memory

1. Client keeps the messages of the current chat.
2. Server sends the earlier messages to the LLM with the new question.
3. Limit how many old messages we send, so the prompt stays small.

Done when: follow-up questions like "and after that?" work.

## Step 3 - Suggested documents

1. Prompt: tell the LLM to write each suggested document as JSON in a fixed form
   (type, name, data format, description).
2. Server: find the JSON in the answer and check it. Drop anything invalid.
3. Client: show each valid suggestion as a card with a "Create" button.
4. Server: on "Create", make the document with normal Hera code.

Done when: the LLM suggests a document, the user clicks, and the document appears in the project.

## Later (not now)

- Generate notebooks (for example, a notebook that runs a simulation).
- Add the project's own documents to the RAG index.
- Save chats per project.

## Open questions

- Which LLM model is good enough?
- Which document types may the LLM suggest? Start with a small list?
- Does the chat replace the right panel, or sit next to it?
